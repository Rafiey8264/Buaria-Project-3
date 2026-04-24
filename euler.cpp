/*
 * ME 5302 Project 3 — 2D Euler Solver (Diamond Airfoil)
 * AUSM flux splitting + MUSCL-TVD reconstruction
 * Rewritten from scratch to match updated project PDF (proj3-3.pdf)
 *
 * Compile: g++ -O2 -o euler euler.cpp -lm
 * Run:     ./euler input.txt
 */
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

// ── Global arrays ──────────────────────────────────────────────────
int NX, NY;
std::vector<std::vector<double>> X, Y;                    // grid coords
std::vector<std::vector<double>> XiX, XiY, EtaX, EtaY, J; // metrics
std::vector<std::vector<std::vector<double>>> Qo, Qn;     // state vectors [4][NX][NY]

// ── Freestream & parameters ────────────────────────────────────────
double rho0, u0, v0, M0, gam, p0, Et0, a0;
double CFL_num;  int maxiter;  double tol;  int muscl_scheme;
std::string gridfile;

// ── Helpers ─────────────────────────────────────────────────────────
static inline double sgn(double x){ return (x>0)-(x<0); }
static inline double minmod(double a, double b){
    if(a*b<=0.0) return 0.0;
    return (std::fabs(a)<std::fabs(b)) ? a : b;
}
// Strip Windows \r
static void strip(std::string &s){
    while(!s.empty() && (s.back()=='\r'||s.back()=='\n')) s.pop_back();
}

// ── READ INPUT ─────────────────────────────────────────────────────
void readInput(const char* fn){
    std::ifstream f(fn); if(!f){ std::cerr<<"Cannot open "<<fn<<"\n"; exit(1);}
    std::string line,key;
    while(std::getline(f,line)){
        strip(line); if(line.empty()||line[0]=='#') continue;
        std::istringstream ss(line); ss>>key;
        if     (key=="rho0")         ss>>rho0;
        else if(key=="u0")           ss>>u0;
        else if(key=="v0")           ss>>v0;
        else if(key=="M0")           ss>>M0;
        else if(key=="gamma")        ss>>gam;
        else if(key=="CFL")          ss>>CFL_num;
        else if(key=="maxiter")      ss>>maxiter;
        else if(key=="tol")          ss>>tol;
        else if(key=="muscl_scheme") ss>>muscl_scheme;
        else if(key=="gridfile")     { ss>>gridfile; strip(gridfile); }
    }
    double Vmag=std::sqrt(u0*u0+v0*v0);
    a0 = Vmag/M0;
    p0 = rho0*a0*a0/gam;
    Et0= p0/(gam-1.0)+0.5*rho0*(u0*u0+v0*v0);
    std::cout<<"Freestream: rho="<<rho0<<" u="<<u0<<" v="<<v0
             <<" M="<<M0<<" p="<<p0<<" Et="<<Et0<<"\n";
}

// ── READ GRID ──────────────────────────────────────────────────────
void readGrid(){
    std::ifstream f(gridfile);
    if(!f){ std::cerr<<"Cannot open grid ["<<gridfile<<"]\n"; exit(1);}
    std::string line; std::getline(f,line); strip(line);
    for(auto &c:line) if(c==',') c=' ';
    std::istringstream(line)>>NX>>NY;
    std::cout<<"Grid: "<<NX<<" x "<<NY<<"\n";
    X.assign(NX,std::vector<double>(NY));
    Y.assign(NX,std::vector<double>(NY));
    for(int i=0;i<NX;i++){
        // skip blank lines until we find data
        while(std::getline(f,line)){
            strip(line);
            bool hasData=false;
            for(char c:line) if(!std::isspace(c)){hasData=true;break;}
            if(hasData){
                std::istringstream ls(line); ls>>X[i][0]>>Y[i][0];
                break;
            }
        }
        for(int j=1;j<NY;j++) f>>X[i][j]>>Y[i][j];
        std::getline(f,line); // consume newline
    }
    std::cout<<"x:["<<X[0][0]<<","<<X[NX-1][0]<<"] y:["<<Y[0][0]<<","<<Y[0][NY-1]<<"]\n";
}

// ── GRID METRICS (CD2 interior, FD2/BD2 boundary) ──────────────────
void gridMetrics(){
    XiX.assign(NX,std::vector<double>(NY));
    XiY.assign(NX,std::vector<double>(NY));
    EtaX.assign(NX,std::vector<double>(NY));
    EtaY.assign(NX,std::vector<double>(NY));
    J.assign(NX,std::vector<double>(NY));
    for(int i=0;i<NX;i++) for(int j=0;j<NY;j++){
        double xxi,yxi,xet,yet;
        // xi derivatives
        if(i==0){
            xxi=(-3*X[0][j]+4*X[1][j]-X[2][j])/2.0;
            yxi=(-3*Y[0][j]+4*Y[1][j]-Y[2][j])/2.0;
        }else if(i==NX-1){
            xxi=(3*X[i][j]-4*X[i-1][j]+X[i-2][j])/2.0;
            yxi=(3*Y[i][j]-4*Y[i-1][j]+Y[i-2][j])/2.0;
        }else{
            xxi=(X[i+1][j]-X[i-1][j])/2.0;
            yxi=(Y[i+1][j]-Y[i-1][j])/2.0;
        }
        // eta derivatives
        if(j==0){
            xet=(-3*X[i][0]+4*X[i][1]-X[i][2])/2.0;
            yet=(-3*Y[i][0]+4*Y[i][1]-Y[i][2])/2.0;
        }else if(j==NY-1){
            xet=(3*X[i][j]-4*X[i][j-1]+X[i][j-2])/2.0;
            yet=(3*Y[i][j]-4*Y[i][j-1]+Y[i][j-2])/2.0;
        }else{
            xet=(X[i][j+1]-X[i][j-1])/2.0;
            yet=(Y[i][j+1]-Y[i][j-1])/2.0;
        }
        double det = xxi*yet - xet*yxi;
        J[i][j]    = 1.0/det;
        XiX[i][j]  =  J[i][j]*yet;
        XiY[i][j]  = -J[i][j]*xet;
        EtaX[i][j] = -J[i][j]*yxi;
        EtaY[i][j] =  J[i][j]*xxi;
    }
}

// ── INIT ───────────────────────────────────────────────────────────
void init(){
    Qo.assign(4,std::vector<std::vector<double>>(NX,std::vector<double>(NY)));
    Qn=Qo;
    for(int i=0;i<NX;i++) for(int j=0;j<NY;j++){
        Qo[0][i][j]=rho0; Qo[1][i][j]=rho0*u0;
        Qo[2][i][j]=rho0*v0; Qo[3][i][j]=Et0;
    }
    Qn=Qo;
}

// ── PRIMITIVES from Q ──────────────────────────────────────────────
struct Prim{ double rho,u,v,p,Et; };
inline Prim prim(const std::vector<std::vector<std::vector<double>>>& Q, int i, int j){
    Prim w;
    w.rho=Q[0][i][j]; w.u=Q[1][i][j]/w.rho; w.v=Q[2][i][j]/w.rho;
    w.Et=Q[3][i][j]; w.p=(gam-1.0)*(w.Et-0.5*w.rho*(w.u*w.u+w.v*w.v));
    return w;
}

// ── BOUNDARY CONDITIONS ────────────────────────────────────────────
// Inlet: i=0, freestream
void inletBC(std::vector<std::vector<std::vector<double>>>& Q){
    for(int j=0;j<NY;j++){
        Q[0][0][j]=rho0; Q[1][0][j]=rho0*u0;
        Q[2][0][j]=rho0*v0; Q[3][0][j]=Et0;
    }
}

// Outlet: i=NX-1, extrapolate
void outletBC(std::vector<std::vector<std::vector<double>>>& Q){
    for(int j=0;j<NY;j++) for(int k=0;k<4;k++)
        Q[k][NX-1][j]=Q[k][NX-2][j];
}

// Wall: j=0 and j=NY-1, free-slip
// Implements all 4 conditions from PDF Section 8
void wallBC(std::vector<std::vector<std::vector<double>>>& Q){
    for(int wall=0;wall<2;wall++){
        int jw = (wall==0) ? 0 : NY-1;
        int ji = (wall==0) ? 1 : NY-2;
        for(int i=0;i<NX;i++){
            // Interior primitives
            Prim wi = prim(Q, i, ji);

            // Wall metrics at wall point
            double sx = XiX[i][jw],  sy = XiY[i][jw];    // xi direction (tangent)
            double ex = EtaX[i][jw], ey = EtaY[i][jw];   // eta direction (normal)
            double smag = std::sqrt(sx*sx + sy*sy);
            double emag = std::sqrt(ex*ex + ey*ey);

            // Normalized metric components (PDF Eq 10)
            double sxh = sx / smag;
            double syh = sy / smag;
            double exh = ex / emag;
            double eyh = ey / emag;

            // PDF Condition 2: d(U_tilde)/d_eta = 0
            // → wall's normalized contravariant tangent velocity = interior's
            double Uti = sxh*wi.u + syh*wi.v;

            // PDF Condition 1: V_tilde = 0 (no penetration)
            // Solve 2x2 system for (uw, vw):
            //   sxh*uw + syh*vw = Uti        (tangential = interior)
            //   exh*uw + eyh*vw = 0          (no penetration)
            double det = sxh*eyh - syh*exh;

            double uw, vw;
            if(std::fabs(det) < 1e-14){
                // Metric determinant degenerate — fall back to projection
                double Vn = wi.u*exh + wi.v*eyh;
                uw = wi.u - Vn*exh;
                vw = wi.v - Vn*eyh;
            } else {
                // Cramer's rule
                uw = ( Uti*eyh - syh*0.0 ) / det;
                vw = ( sxh*0.0 - Uti*exh ) / det;
            }

            // PDF Condition 3: dP/d_eta = 0 at wall → p_wall = p_interior
            double pw = wi.p;

            // PDF Condition 4: stagnation enthalpy conserved
            double h0_int = gam/(gam-1.0)*wi.p/wi.rho + 0.5*(wi.u*wi.u + wi.v*wi.v);
            double KE_w = 0.5*(uw*uw + vw*vw);
            double denom = h0_int - KE_w;

            double rhow;
            if(denom <= 0.0){
                rhow = wi.rho;  // fallback
            } else {
                rhow = gam * pw / ((gam-1.0) * denom);
            }

            double Etw = pw/(gam-1.0) + 0.5*rhow*(uw*uw+vw*vw);

            Q[0][i][jw]=rhow;
            Q[1][i][jw]=rhow*uw;
            Q[2][i][jw]=rhow*vw;
            Q[3][i][jw]=Etw;
        }
    }
}

// ── MUSCL reconstruction (updated PDF Section 5.C) ─────────────────
// Returns QL, QR at interface ic+1/2 (dir=0: xi, dir=1: eta)
void muscl(int i, int j, int dir, double QL[4], double QR[4],
           const std::vector<std::vector<std::vector<double>>>& Q){

    // Stencil: im1, ic, ip1, ip2
    int im1_i=i, im1_j=j, ic_i=i, ic_j=j;
    int ip1_i=i, ip1_j=j, ip2_i=i, ip2_j=j;
    bool at_boundary = false;

    if(dir==0){
        im1_i=i-1; ic_i=i; ip1_i=i+1; ip2_i=i+2;
        if(i-1<0 || i+2>NX-1) at_boundary=true;
    } else {
        im1_j=j-1; ic_j=j; ip1_j=j+1; ip2_j=j+2;
        if(j-1<0 || j+2>NY-1) at_boundary=true;
    }

    // Clamp
    auto cl=[](int v,int lo,int hi){return std::max(lo,std::min(hi,v));};
    im1_i=cl(im1_i,0,NX-1); im1_j=cl(im1_j,0,NY-1);
    ic_i =cl(ic_i, 0,NX-1); ic_j =cl(ic_j, 0,NY-1);
    ip1_i=cl(ip1_i,0,NX-1); ip1_j=cl(ip1_j,0,NY-1);
    ip2_i=cl(ip2_i,0,NX-1); ip2_j=cl(ip2_j,0,NY-1);

    for(int k=0;k<4;k++){
        double qi  = Q[k][ic_i][ic_j];
        double qp1 = Q[k][ip1_i][ip1_j];

        if(muscl_scheme==0 || (muscl_scheme==2 && at_boundary)){
            // Case (a): first order; also case (c) at boundaries
            QL[k]=qi; QR[k]=qp1;
        } else if(muscl_scheme==1){
            // Case (b): PURE central average as stated in PDF Eq 21
            // Q_L = Q_R = (Q_i + Q_{i+1}) / 2
            // Note: this scheme may produce limit-cycle or divergence due to
            // Godunov's theorem (linear non-upwind scheme cannot be monotone).
            QL[k] = 0.5*(qi + qp1);
            QR[k] = 0.5*(qi + qp1);
        } else if(muscl_scheme==3){
            // Case (b) STABILIZED: central average with adaptive dissipation.
            // When local variation exceeds a small threshold (shock sensor),
            // fall back to 1st-order upwind to stabilize Scheme 1.
            // This is NOT the PDF Scheme 1 — it is a stabilized variant.
            double Qim1 = (at_boundary) ? qi  : Q[k][im1_i][im1_j];

            bool use_upwind = false;
            if(!at_boundary){
                double Qmax = std::max({std::fabs(Qim1), std::fabs(qi),
                                       std::fabs(qp1)}) + 1e-12;
                double var = std::fabs(qp1 - qi) / Qmax;
                if(var > 0.0001) use_upwind = true;
            }

            if(use_upwind){
                QL[k] = qi;
                QR[k] = qp1;
            } else {
                QL[k] = 0.5*(qi + qp1);
                QR[k] = 0.5*(qi + qp1);
            }
        } else {
            // Case (c): TVD with Chakravarthy-Osher, Eq 22-26
            double qm1 = Q[k][im1_i][im1_j];
            double qp2 = Q[k][ip2_i][ip2_j];
            double dm = qi  - qm1;   // ΔQ_{i-1/2}
            double dp = qp1 - qi;    // ΔQ_{i+1/2}
            double d3 = qp2 - qp1;   // ΔQ_{i+3/2}

            double phiL=0.0, phiR=0.0;
            if(std::fabs(dm)>1e-30) phiL = minmod(dm, 2.0*dp)/dm;
            if(std::fabs(d3)>1e-30) phiR = minmod(d3, 2.0*dp)/d3;

            QL[k] = qi  + 0.5*dm*phiL;   // Eq 22 (plus sign, confirmed)
            QR[k] = qp1 - 0.5*d3*phiR;   // Eq 23
        }
    }
}

// ── AUSM FLUX in ξ direction: F'_{i+1/2,j} ─────────────────────────
void fluxFX(int i, int j, double fx[4],
            const std::vector<std::vector<std::vector<double>>>& Q){
    double QL[4], QR[4];
    muscl(i, j, 0, QL, QR, Q);

    // L,R primitives (with positivity check)
    double rL=QL[0], uL=QL[1]/rL, vL=QL[2]/rL, EtL=QL[3];
    double pL=(gam-1)*(EtL-0.5*rL*(uL*uL+vL*vL));
    double rR=QR[0], uR=QR[1]/rR, vR=QR[2]/rR, EtR=QR[3];
    double pR=(gam-1)*(EtR-0.5*rR*(uR*uR+vR*vR));

    // Warn on non-physical states but continue (many CFD codes clamp)
    static int warn_count = 0;
    if((rL<=0 || pL<=0 || rR<=0 || pR<=0) && warn_count < 10){
        std::cerr << "WARNING: non-physical state in flux at i=" << i
                  << " j=" << j << " pL=" << pL << " pR=" << pR
                  << " rL=" << rL << " rR=" << rR << "\n";
        warn_count++;
    }

    // Clamp to positive values to avoid NaN in sqrt
    double rL_safe = std::max(1e-10, rL);
    double pL_safe = std::max(1e-10, pL);
    double rR_safe = std::max(1e-10, rR);
    double pR_safe = std::max(1e-10, pR);
    double aL = std::sqrt(gam*pL_safe/rL_safe);
    double aR = std::sqrt(gam*pR_safe/rR_safe);

    // Averaged metrics at i+1/2
    double sxb=0.5*(XiX[i][j]+XiX[i+1][j]);
    double syb=0.5*(XiY[i][j]+XiY[i+1][j]);
    double Jb =0.5*(J[i][j]+J[i+1][j]);
    double smag=std::sqrt(sxb*sxb+syb*syb);

    // Normalized contravariant velocity (Eq 10)
    double UtL=(sxb*uL+syb*vL)/smag;
    double UtR=(sxb*uR+syb*vR)/smag;

    // Mach (Eq 11)
    double ML=UtL/aL, MR=UtR/aR;

    // Split Mach & pressure (Eq 12-15)
    double Mp,Mm,pp,pm;
    if(std::fabs(ML)<=1){ Mp=0.25*(ML+1)*(ML+1); pp=(2-ML)*Mp; }
    else { Mp=0.5*(ML+std::fabs(ML)); pp=0.5*(1+sgn(ML)); }
    if(std::fabs(MR)<=1){ Mm=-0.25*(MR-1)*(MR-1); pm=-(2+MR)*Mm; }
    else { Mm=0.5*(MR-std::fabs(MR)); pm=0.5*(1-sgn(MR)); }

    // Contravariant velocity contributions (Eq 16)
    double Up=Mp*aL*smag, Um=Mm*aR*smag;

    // Flux (Eq 17) — use safe pressures for consistency with sound speed
    double Jinv=1.0/Jb;
    fx[0]=Jinv*(rL_safe*Up + rR_safe*Um);
    fx[1]=Jinv*(rL_safe*uL*Up + rR_safe*uR*Um + sxb*(pp*pL_safe+pm*pR_safe));
    fx[2]=Jinv*(rL_safe*vL*Up + rR_safe*vR*Um + syb*(pp*pL_safe+pm*pR_safe));
    fx[3]=Jinv*(Up*(EtL+pL_safe) + Um*(EtR+pR_safe));
}

// ── AUSM FLUX in η direction: G'_{i,j+1/2} ─────────────────────────
void fluxGY(int i, int j, double gy[4],
            const std::vector<std::vector<std::vector<double>>>& Q){
    double QL[4], QR[4];
    muscl(i, j, 1, QL, QR, Q);

    double rL=QL[0], uL=QL[1]/rL, vL=QL[2]/rL, EtL=QL[3];
    double pL=(gam-1)*(EtL-0.5*rL*(uL*uL+vL*vL));
    double rR=QR[0], uR=QR[1]/rR, vR=QR[2]/rR, EtR=QR[3];
    double pR=(gam-1)*(EtR-0.5*rR*(uR*uR+vR*vR));

    static int warn_count_g = 0;
    if((rL<=0 || pL<=0 || rR<=0 || pR<=0) && warn_count_g < 10){
        std::cerr << "WARNING: non-physical state in fluxGY at i=" << i
                  << " j=" << j << " pL=" << pL << " pR=" << pR
                  << " rL=" << rL << " rR=" << rR << "\n";
        warn_count_g++;
    }

    double rL_safe = std::max(1e-10, rL);
    double pL_safe = std::max(1e-10, pL);
    double rR_safe = std::max(1e-10, rR);
    double pR_safe = std::max(1e-10, pR);
    double aL = std::sqrt(gam*pL_safe/rL_safe);
    double aR = std::sqrt(gam*pR_safe/rR_safe);

    // Averaged metrics at j+1/2
    double exb=0.5*(EtaX[i][j]+EtaX[i][j+1]);
    double eyb=0.5*(EtaY[i][j]+EtaY[i][j+1]);
    double Jb =0.5*(J[i][j]+J[i][j+1]);
    double emag=std::sqrt(exb*exb+eyb*eyb);

    // Normalized V_tilde (Eq 10, eta version)
    double VtL=(exb*uL+eyb*vL)/emag;
    double VtR=(exb*uR+eyb*vR)/emag;

    // Mach from V_tilde (Eq 18)
    double ML=VtL/aL, MR=VtR/aR;

    double Mp,Mm,pp,pm;
    if(std::fabs(ML)<=1){ Mp=0.25*(ML+1)*(ML+1); pp=(2-ML)*Mp; }
    else { Mp=0.5*(ML+std::fabs(ML)); pp=0.5*(1+sgn(ML)); }
    if(std::fabs(MR)<=1){ Mm=-0.25*(MR-1)*(MR-1); pm=-(2+MR)*Mm; }
    else { Mm=0.5*(MR-std::fabs(MR)); pm=0.5*(1-sgn(MR)); }

    // Eq 19
    double Vp=Mp*aL*emag, Vm=Mm*aR*emag;

    // Eq 20 — use safe pressures for consistency with sound speed
    double Jinv=1.0/Jb;
    gy[0]=Jinv*(rL_safe*Vp + rR_safe*Vm);
    gy[1]=Jinv*(rL_safe*uL*Vp + rR_safe*uR*Vm + exb*(pp*pL_safe+pm*pR_safe));
    gy[2]=Jinv*(rL_safe*vL*Vp + rR_safe*vR*Vm + eyb*(pp*pL_safe+pm*pR_safe));
    gy[3]=Jinv*(Vp*(EtL+pL_safe) + Vm*(EtR+pR_safe));
}

// ── CALC DT (Eq 27-28) ────────────────────────────────────────────
double calcDt(){
    double mx=0;
    static int dt_warn = 0;
    for(int i=0;i<NX;i++) for(int j=0;j<NY;j++){
        Prim w=prim(Qo,i,j);

        // Positivity check: warn (don't abort) if non-physical state
        if((w.rho <= 0.0 || w.p <= 0.0) && dt_warn < 5){
            std::cerr << "WARNING: non-physical state in calcDt at i=" << i
                      << " j=" << j << " rho=" << w.rho << " p=" << w.p << "\n";
            dt_warn++;
        }

        // Clamp for A^2 computation only (so we can see where it fails)
        double p_safe = std::max(1e-10, w.p);
        double rho_safe = std::max(1e-10, w.rho);

        double U=XiX[i][j]*w.u+XiY[i][j]*w.v;
        double V=EtaX[i][j]*w.u+EtaY[i][j]*w.v;
        double A2=(gam*p_safe/rho_safe)*(
            (XiX[i][j]*XiX[i][j]+XiY[i][j]*XiY[i][j])+
            (EtaX[i][j]*EtaX[i][j]+EtaY[i][j]*EtaY[i][j])+
            2*std::fabs(XiX[i][j]*EtaX[i][j]+XiY[i][j]*EtaY[i][j]));
        double val=std::fabs(U)+std::fabs(V)+std::sqrt(A2);
        if(val>mx) mx=val;
    }
    return CFL_num/mx;
}

// ── SOLVE (Eq 8) ──────────────────────────────────────────────────
void solve(double dt){
    for(int j=1;j<NY-1;j++) for(int i=1;i<NX-1;i++){
        double FXp[4],FXm[4],GYp[4],GYm[4];
        fluxFX(i,  j,FXp,Qo);   // F'_{i+1/2}
        fluxFX(i-1,j,FXm,Qo);   // F'_{i-1/2}
        fluxGY(i,j,  GYp,Qo);   // G'_{j+1/2}
        fluxGY(i,j-1,GYm,Qo);   // G'_{j-1/2}
        for(int k=0;k<4;k++)
            Qn[k][i][j]=Qo[k][i][j]-dt*J[i][j]*((FXp[k]-FXm[k])+(GYp[k]-GYm[k]));
    }
}

// ── ERROR (Eq 31-32) ──────────────────────────────────────────────
double calcError(){
    double s=0; int Nt=(NX-2)*(NY-2);
    for(int i=1;i<NX-1;i++) for(int j=1;j<NY-1;j++)
        for(int k=0;k<4;k++){ double d=Qn[k][i][j]-Qo[k][i][j]; s+=d*d; }
    return std::sqrt(s/(4.0*Nt));
}

// ── OUTPUT ─────────────────────────────────────────────────────────
void writeOutput(const std::string& pre){
    std::string fn=pre+"_solution.dat";
    std::ofstream f(fn);
    f<<NX<<" "<<NY<<"\n";
    for(int i=0;i<NX;i++) for(int j=0;j<NY;j++){
        Prim w=prim(Qn,i,j);
        double a=std::sqrt(gam*std::max(1e-10, w.p)/std::max(1e-10, w.rho));
        double M=std::sqrt(w.u*w.u+w.v*w.v)/a;
        f<<std::scientific<<std::setprecision(10)
         <<X[i][j]<<" "<<Y[i][j]<<" "<<w.rho<<" "<<w.u<<" "<<w.v
         <<" "<<w.p<<" "<<M<<"\n";
    }
    std::cout<<"Written "<<fn<<"\n";
}
void writeGrid(){
    std::ofstream f("grid_output.dat");
    f<<NX<<" "<<NY<<"\n";
    for(int i=0;i<NX;i++) for(int j=0;j<NY;j++)
        f<<std::scientific<<std::setprecision(10)<<X[i][j]<<" "<<Y[i][j]<<"\n";
}
void writeConv(const std::vector<double>& err, const std::string& pre){
    std::ofstream f(pre+"_convergence.dat");
    for(size_t n=0;n<err.size();n++)
        f<<n+1<<" "<<std::scientific<<std::setprecision(10)<<err[n]<<"\n";
}

// ── MAIN ──────────────────────────────────────────────────────────
int main(int argc, char* argv[]){
    if(argc<2){ std::cerr<<"Usage: ./euler input.txt\n"; return 1; }

    readInput(argv[1]);
    readGrid();
    gridMetrics();
    init();
    writeGrid();

    double time=0, E1=0;
    std::vector<double> rel_err;

    std::cout<<std::setw(8)<<"Step"<<std::setw(16)<<"Time"
             <<std::setw(16)<<"dt"<<std::setw(16)<<"AbsErr"
             <<std::setw(16)<<"RelErr"<<"\n";

    for(int n=1;n<=maxiter;n++){
        // ── Code outline order ──
        // BCs on Qold before solve
        inletBC(Qo);
        wallBC(Qo);

        // CALCDT
        double dt=calcDt();
        time+=dt;

        // SOLVE: Qn interior from Qo
        solve(dt);

        // BCs on Qnew after solve (so Qnew is complete before Qo=Qn)
        inletBC(Qn);
        wallBC(Qn);
        outletBC(Qn);

        // QERROR
        double E=calcError();
        if(n==1) E1=E;
        double R=(E1>0)?E/E1:0;
        rel_err.push_back(R);

        if(n<=10||n%500==0||R<tol)
            std::cout<<std::setw(8)<<n
                     <<std::setw(16)<<std::scientific<<std::setprecision(6)<<time
                     <<std::setw(16)<<dt<<std::setw(16)<<E<<std::setw(16)<<R<<"\n";

        if(R<tol&&n>1){
            std::cout<<"\nConverged at step "<<n<<" R="<<R<<"\n";
            break;
        }

        // Qold = Qnew
        Qo=Qn;
    }

    std::string pre="scheme"+std::to_string(muscl_scheme);
    writeOutput(pre);
    writeConv(rel_err,pre);
    std::cout<<"Done!\n";
    return 0;
}

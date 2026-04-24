
The grid files can be read in following manner:

----------------

  OPEN file 

  READ:  nx, ny

  ALLOCATE arrays  x(nx, ny), y(nx, ny)

  loopi: for i = 1 : nx

    READ:  <blank>

    loopj: for j = 1 : ny
      READ:  x(i,j), y(i,j) 
    end loopj

  end loopi

  CLOSE file

----------------


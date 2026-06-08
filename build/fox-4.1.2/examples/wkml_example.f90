program wkml_example
! this program read car accidents data of Cambridge 

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  integer i,m
  parameter(m=10) 

  double precision :: latitude(m), longitude(m)  


  open(unit=90, file='wkml_example_input.txt',status='old')


10 FORMAT(2x,F9.6,5x,F8.6)
  do i=1,m
    read (90,10) latitude(i), longitude(i)
  end do

!   print*, longitude

  close(90)


  call kmlBeginFile(myfile, "wkml_example.kml", -1)
  call kmlCreatePoints(myfile, longitude, latitude)
  call kmlFinishFile(myfile)


end program wkml_example

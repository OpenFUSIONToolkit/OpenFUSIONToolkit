program test_sax_reader

  use m_sax_reader

  type(file_buffer_t) :: fb

  integer :: iostat
  character(1) :: c

  call open_file(fb, file="test_sax_reader_1.in",  iostat=iostat)

  c = read_char(fb, iostat)
  write(*,'(2a)') 'char:', c
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader

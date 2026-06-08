program test_sax_reader

  use m_sax_reader

  type(file_buffer_t) :: fb

  character(len=1) :: c
  character(len=3) :: s
  integer :: iostat

  call open_file(fb, file="test_sax_reader_3.in", iostat=iostat)

  s = read_chars(fb, 3, iostat)

  call push_chars(fb, 'x')

  c = read_char(fb, iostat)

  write(*,'(3a)') 'char:', s, ':'
  write(*,'(3a)') 'char:', c, ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader

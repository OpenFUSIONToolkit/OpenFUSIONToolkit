program test_sax_reader

  use m_sax_reader

  type(file_buffer_t) :: fb

  character(len=5) :: s
  integer :: iostat

  call open_file(fb, file="test_sax_reader_5.in", iostat=iostat)

  call push_chars(fb, "abc")

  s = read_chars(fb, 5, iostat)
  write(*,'(3a)') 'char:', s, ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader

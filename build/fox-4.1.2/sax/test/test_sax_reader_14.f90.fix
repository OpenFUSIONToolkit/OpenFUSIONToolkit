program test_sax_reader

  use m_sax_reader

  character :: c
  type(file_buffer_t) :: fb

  integer :: iostat

! test 1: quick & easy

  call open_file(fb, file="test_sax_reader_14.in", iostat=iostat)

  call get_characters_until_all_of(fb, '--', iostat) 

  write(*,'(3a)') 'char:', retrieve_namebuffer(fb), ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

! test two - string over 1000 chars

  call open_file(fb, file="test_sax_reader_14b.in", iostat=iostat) 

  call get_characters_until_all_of(fb, '--', iostat) 

  write(*,'(3a)') 'char:', retrieve_namebuffer(fb), ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader

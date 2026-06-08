program test_sax_reader

  use m_sax_reader

  character(40) :: buf
  character :: c
  type(file_buffer_t) :: fb

  integer :: iostat

  ! test 1: move 40 chars in, and then pick up the name

  call open_file(fb, file="test_sax_reader_12a.in", iostat=iostat)

  buf = get_characters(fb, 40, iostat) 

  call get_characters_until_not_one_of(fb, 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ', iostat) 

  write(*,'(3a)') 'char:', retrieve_namebuffer(fb), ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

  ! test 2: move 40 chars in, and then pick up the name, which is now 2004 chars long.

  call open_file(fb, file="test_sax_reader_12b.in", iostat=iostat) 

  buf = get_characters(fb, 40, iostat) 

  call get_characters_until_not_one_of(fb, 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ', iostat) 

  print*, len_namebuffer(fb)

  write(*,'(3a)') 'char:', retrieve_namebuffer(fb), ':'
  write(*,'(a,i0)') 'iost:', iostat

  call get_characters_until_not_one_of(fb, 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ', iostat) 

  print*, len_namebuffer(fb)

  write(*,'(3a)') 'char:', retrieve_namebuffer(fb), ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader

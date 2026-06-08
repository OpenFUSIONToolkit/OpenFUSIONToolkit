module checksub
contains
function check(c) result(p)
character, intent(in) :: c
logical :: p
p = (c=='c')
end function check
end module checksub

program test_sax_reader

  use checksub
  use m_sax_reader

  character :: c
  type(file_buffer_t) :: fb

  integer :: iostat

! test 1: quick & easy

  call open_file(fb, string="abcdef", iostat=iostat)

  call get_characters_until_condition(fb, check, .true., iostat) 

  write(*,'(3a)') 'char:', retrieve_namebuffer(fb), ':'
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader

program test_sax_reader

  use m_common_array_str
  use m_sax_parser
  use m_sax_reader
  use m_sax_types

  type(file_buffer_t) :: fb
  type(sax_parser_t) :: fx

  integer :: iostat

  call open_file(fb, file="test_sax_fsm_1.in", iostat=iostat)

  call sax_parser_init(fx, "test_sax_fsm_1.in")

  call sax_parse(fx, fb)

  write(*,'(i0)') fx%xds%xml_version
  write(*,'(a)') str_vs(fx%xds%encoding)
  write(*,'(l1)') fx%xds%standalone

  call sax_parser_destroy(fx)
  call close_file(fb) 

end program test_sax_reader

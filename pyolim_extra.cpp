// Extra definitions that aren't generated

py::enum_<state>(m, "State")
  .value("Valid", state::valid)
  .value("Trial", state::trial)
  .value("Far", state::far);

py::class_<fac_src>(m, "FacCenter")
  .def(py::init<double, double, double>())
  .def_readwrite("i", &fac_src::i)
  .def_readwrite("j", &fac_src::j)
  .def_readwrite("s", &fac_src::s)
  ;

py::class_<fac_src_3d>(m, "FacCenter3d")
  .def(py::init<double, double, double, double>())
  .def_readwrite("i", &fac_src_3d::i)
  .def_readwrite("j", &fac_src_3d::j)
  .def_readwrite("k", &fac_src_3d::k)
  .def_readwrite("s", &fac_src_3d::s)
  ;

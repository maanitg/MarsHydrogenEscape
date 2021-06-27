#include "MarsHydrogenEscape.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

PYBIND11_MODULE(Mars, m) {
  py::class_<MarsHydrogenEscape>(m, "MarsHydrogenEscape")
    .def(py::init())
    .def("WarmTimeAfterSeveralImpacts", 
         &MarsHydrogenEscape::WarmTimeAfterSeveralImpacts)
    .def_property("S_1", 
           &MarsHydrogenEscape::get_S_1, 
           &MarsHydrogenEscape::set_S_1)
    .def_property("grav", 
           &MarsHydrogenEscape::get_grav, 
           &MarsHydrogenEscape::set_grav)
    .def_property("area", 
           &MarsHydrogenEscape::get_area, 
           &MarsHydrogenEscape::set_area)
    .def_property("rtol", 
           &MarsHydrogenEscape::get_rtol, 
           &MarsHydrogenEscape::set_rtol)
    .def_property("atol", 
           &MarsHydrogenEscape::get_atol, 
           &MarsHydrogenEscape::set_atol)
    .def_property("impact_diameters", 
           [](MarsHydrogenEscape &self)
           {
             py::array out = py::cast(self.get_impact_diameters());
             return out;
           },
           &MarsHydrogenEscape::set_impact_diameters);
  
}
#include "functions.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

PYBIND11_MODULE(dispersion, m) {
    m.doc() = "Module for dispersion curves generation for P-SV elastic guided waves";

    m.def("compound_Q_matrix",
    [](py::array_t<std::complex<data_t>, py::array::c_style> Q, data_t mu, std::complex<data_t> r, std::complex<data_t> s, data_t t)
    {compound_Q_matrix(Q.mutable_data(),mu,r,s,t);},
    "Build the compound Q matrix.",
    py::arg("matrix"), py::arg("mu"), py::arg("r"), py::arg("s"), py::arg("t"));

    m.def("compound_invQ_matrix",
    [](py::array_t<std::complex<data_t>, py::array::c_style> Q, data_t mu, std::complex<data_t> r, std::complex<data_t> s, data_t g)
    {compound_invQ_matrix(Q.mutable_data(),mu,r,s,g);},
    "Build the compound inverse Q matrix.",
    py::arg("matrix"), py::arg("mu"), py::arg("r"), py::arg("s"), py::arg("g"));

    m.def("compound_field_matrix",
    [](py::array_t<std::complex<data_t>, py::array::c_style> Q, data_t f, data_t k, data_t a, data_t b, data_t rho, data_t z)
    {compound_field_matrix(Q.mutable_data(),f,k,a,b,rho,z);},
    "Build the compound field matrix.",
    py::arg("matrix"), py::arg("frequency"), py::arg("wavenumber"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("compound_inv_field_matrix",
    [](py::array_t<std::complex<data_t>, py::array::c_style> Q, data_t f, data_t k, data_t a, data_t b, data_t rho, data_t z)
    {compound_inv_field_matrix(Q.mutable_data(),f,k,a,b,rho,z);},
    "Build the compound inverse field matrix.",
    py::arg("matrix"), py::arg("frequency"), py::arg("wavenumber"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("compound_propagator_matrix",
    [](py::array_t<std::complex<data_t>, py::array::c_style> Q, data_t f, data_t k, data_t a, data_t b, data_t rho, data_t z)
    {compound_propagator_matrix(Q.mutable_data(),f,k,a,b,rho,z);},
    "Build the compound propagator matrix.",
    py::arg("matrix"), py::arg("frequency"), py::arg("wavenumber"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("free_solid_free",
    [](const data_t f, const data_t k, const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d) -> std::complex<data_t>
    {return free_solid_free(f,k,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for Lamb waves (elastic plate in a free space).",
    py::arg("frequency"), py::arg("wavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("free_solid_free_map",
    [](py::array_t<std::complex<data_t>, py::array::c_style> val, 
                py::array_t<data_t, py::array::c_style> f,
                const int nf,
                py::array_t<data_t, py::array::c_style> k,
                const int nk, 
                const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d)
    {free_solid_free_map(val.mutable_data(),f.mutable_data(),nf,k.mutable_data(),nk,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for Lamb waves (elastic plate in a free space).",
    py::arg("output"), py::arg("frequency"), py::arg("nfrequency"), py::arg("wavenumber"), py::arg("nwavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("rigid_solid_rigid",
    [](const data_t f, const data_t k, const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d) -> std::complex<data_t>
    {return rigid_solid_rigid(f,k,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for rigid guided waves (elastic layer with rigid walls).",
    py::arg("frequency"), py::arg("wavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));
    
    m.def("rigid_solid_rigid_map",
    [](py::array_t<std::complex<data_t>, py::array::c_style> val, 
                py::array_t<data_t, py::array::c_style> f,
                const int nf,
                py::array_t<data_t, py::array::c_style> k,
                const int nk, 
                const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d)
    {rigid_solid_rigid_map(val.mutable_data(),f.mutable_data(),nf,k.mutable_data(),nk,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for rigid guided waves (elastic layer with rigid walls).",
    py::arg("output"), py::arg("frequency"), py::arg("nfrequency"), py::arg("wavenumber"), py::arg("nwavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("free_solid_halfspace",
    [](const data_t f, const data_t k, const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d) -> std::complex<data_t>
    {return free_solid_halfspace(f,k,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for surface waves.",
    py::arg("frequency"), py::arg("wavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("free_solid_halfspace_map",
    [](py::array_t<std::complex<data_t>, py::array::c_style> val, 
                py::array_t<data_t, py::array::c_style> f,
                const int nf,
                py::array_t<data_t, py::array::c_style> k,
                const int nk, 
                const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d)
    {free_solid_halfspace_map(val.mutable_data(),f.mutable_data(),nf,k.mutable_data(),nk,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for surface waves.",
    py::arg("output"), py::arg("frequency"), py::arg("nfrequency"), py::arg("wavenumber"), py::arg("nwavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("halfspace_solid_halfspace",
    [](const data_t f, const data_t k, const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d) -> std::complex<data_t>
    {return halfspace_solid_halfspace(f,k,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for normal guided waves (layer embedded between two halfspaces).",
    py::arg("frequency"), py::arg("wavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("halfspace_solid_halfspace_map",
    [](py::array_t<std::complex<data_t>, py::array::c_style> val, 
                py::array_t<data_t, py::array::c_style> f,
                const int nf,
                py::array_t<data_t, py::array::c_style> k,
                const int nk, 
                const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d)
    {halfspace_solid_halfspace_map(val.mutable_data(),f.mutable_data(),nf,k.mutable_data(),nk,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for normal guided waves (layer embedded between two halfspaces).",
    py::arg("output"), py::arg("frequency"), py::arg("nfrequency"), py::arg("wavenumber"), py::arg("nwavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));


    m.def("halfspace_solid_halfspace_leaky",
    [](const data_t f, const std::complex<data_t> k, const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d) -> std::complex<data_t>
    {return halfspace_solid_halfspace_leaky(f,k,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for leaky guided waves (layer embedded between two halfspaces). Normal modes can be obtained when Im(k)=0.",
    py::arg("frequency"), py::arg("wavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

    m.def("free_solid_halfspace_leaky",
    [](const data_t f, const std::complex<data_t> k, const int n,
                py::array_t<data_t, py::array::c_style> a,
                py::array_t<data_t, py::array::c_style> b,
                py::array_t<data_t, py::array::c_style> rho,
                py::array_t<data_t, py::array::c_style> d) -> std::complex<data_t>
    {return free_solid_halfspace_leaky(f,k,n,a.mutable_data(),b.mutable_data(),rho.mutable_data(),d.mutable_data());},
    "Characteristic function for leaky surface waves. Normal modes can be obtained when Im(k)=0.",
    py::arg("frequency"), py::arg("wavenumber"), py::arg("nlayers"), py::arg("vp"), py::arg("vs"), py::arg("density"), py::arg("thickness"));

}

    
/* PYBIND11_MODULE(dispersion, m) {
    m.doc() = "Module for dispersion curves generation for P-SV elastic guided waves";
    m.def("Q_matrix",&Q_matrix,"Build the Q matrix.",
    py::arg("matrix"), py::arg("mu"), py::arg("r"), py::arg("s"), py::arg("t"));
} */
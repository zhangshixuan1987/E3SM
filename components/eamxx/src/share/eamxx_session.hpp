#ifndef SCREAM_SESSION_HPP
#define SCREAM_SESSION_HPP

namespace scream {

int get_default_fpes ();
void initialize_eamxx_session(bool print_config = true);
void initialize_eamxx_session (int argc, char **argv, bool print_config = true);

// A version callable from Fortran, which can help
// in case of errors to correctly shut down Kokkos
extern "C" {
void finalize_eamxx_session();
} // extern "C"
} // namespace scream

#endif // SCREAM_SESSION_HPP

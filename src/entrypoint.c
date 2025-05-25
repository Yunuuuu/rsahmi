// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_rsahmi_extendr(void *dll);

void R_init_rsahmi(void *dll) {
    R_init_rsahmi_extendr(dll);
}

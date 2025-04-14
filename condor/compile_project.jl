using Pkg
Pkg.instantiate()
Pkg.update()
Pkg.precompile()

@info "PACKAGE STATUSES:"
@info Pkg.status()

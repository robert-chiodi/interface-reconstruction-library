# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install irl
#
# You can edit this file again by typing:
#
#     spack edit irl
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack.package import *


class Irl(CMakePackage):
    """Interface Reconstruction Library with computational geometry algorithms used in 
       geometric VOF and ALE remapping schemes, among others.    
    """

    homepage = "https://github.com/robert-chiodi/interface-reconstruction-library"
    git      = "https://github.com/robert-chiodi/interface-reconstruction-library.git"

    # Versions
    version('develop', branch='master', submodules=False, preferred=True)

    # Github accounts to notify of package changes
    maintainers = ['robert-chiodi']

    # Variants
    variant('unit', default=False,
            description='Enable Unit Tests')    
    variant('absl', default=True,
            description='Use ABSL for some underlying classes')            
    variant('fortran', default=False,
            description='Build Fortran library interface')

    # Dependencies
    depends_on('cmake', type='build')    
    depends_on('eigen')

    def cmake_args(self):
        return [
            self.define_from_variant('BUILD_TESTING', 'unit'),
            self.define_from_variant('IRL_USE_ABSL', 'absl'),
            self.define_from_variant('IRL_BUILD_FORTRAN', 'fortran')
        ]

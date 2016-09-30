# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import numpy as np

import moldesign as mdt
import pyccc
from moldesign import units as u, compute
from moldesign.interfaces.ambertools import IMAGE
from .base import QMBase


def exports(o):
    __all__.append(o.__name__)
    return o

__all__ = []


class SQMPotential(QMBase):
    DEFAULT_PROPERTIES = ['potential_energy',
                          'wfn',
                          'mulliken']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    THEORIES = ('MNDO MNDO/d AM1 AM1/d PM3 PDDG PDDG/MNDO PDDG/PM3 RM1 '
               'PM3CARB1 PM3-MAIS PM6 DFTB').split()
    FORCE_UNITS = u.hartree / u.bohr

    def __init__(self, **kwargs):
        super(SQMPotential, self).__init__(**kwargs)

    @mdt.utils.kwargs_from(mdt.compute.run_job)
    def calculate(self, requests=None, guess=None, **kwargs):
        inputfile = ['SQM input for %s' % self.mol.name,
                     ' &qmmm',
                     " qm_theory='%s', qmcharge=%d, printcharges=1, maxcyc=0" % (
                         self.params.theory, self.get_formal_charge()),
                     ' /']
        for atom in self.mol.atoms:
            inputfile.append(' {atom.atnum}  {atom.name}  {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}'.format(
                atom=atom, pos=atom.position.value_in(u.angstrom)))

        image = compute.get_image_path(IMAGE)
        job = pyccc.Job(image=image,
                        command='sqm -i mol.in -o mol.out',
                        inputs={'mol.in': '\n'.join(inputfile)},
                        name="sqm single point, %s" % self.mol.name,
                        when_finished=self._parse_results)

        return compute.run_job(job, _return_result=True, **kwargs)

    def _parse_results(self, job):
        result = {}
        lines = iter(job.get_output('mol.out'))
        while True:
            try:
                line = lines.next()
            except StopIteration:
                break

            fields = line.split()
            if fields[:2] == ['QM', 'DIPOLE']:
                # TODO: CHECK UNITS
                result['dipole'] = np.array(map(float,fields[2:5])) * u.q_e * u.angstrom

            if fields == 'Atomic Charges for Step 1'.split():
                assert lines.next().split() == 'Atom Element Mulliken Charge'.split()
                fields = lines.next().split()
                charges = []
                while fields[:4] != 'Total Mulliken Charge ='.split():
                    charges.append(float(fields[2]))
                result['charges'] = charges * u.q_e

            if fields == 'Final MO eigenvalues (eV):'.split():
                mo_energies = []
                fields = lines.next().split()
                while fields:
                    mo_energies.extend(map(float, fields))
                    fields = lines.next().split()
                result['mo_energies'] = mo_energies * u.eV

            if fields[:3] == 'Heat of formation'.split():
                result['potential_energy'] = float(fields[4]) * u.kcalpermol

        return result


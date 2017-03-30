"""
Add trident version to test results.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
from trident.testing import \
    generate_results, \
    test_results_dir

def test_result_generate():
    if not generate_results:
        return
    try:
        import hglib
    except ImportError:
        raise RuntimeError(
            "Generating test results requires python-hglib. " + \
            "Install with: pip install python-hglib")

    repo = hglib.open("..")
    my_hash = repo.identify().decode("utf").strip().split()[0]
    fh = open(os.path.join(test_results_dir, "TRIDENT_VERSION"), "w")
    fh.write("%s\n" % my_hash)
    fh.close()

#!/usr/bin/env python3
import os
import unittest

class TestHometools(unittest.TestCase):
    def test_cli(self):
        from subprocess import Popen, PIPE
        p = Popen('syri -h'.split(), stdout=PIPE, stderr=PIPE)
        out = p.communicate()
        assert 'error' not in out[1].decode()

    # def test_coords_input(self):
    #     from subprocess import Popen, PIPE
    #     from glob import glob
    #     coords = "out.filtered.coords"
    #     delta = "out.filtered.delta"
    #     ref = "tair10_chr4.fa.gz"
    #     qry = "ler_chr4.fa.gz"
    #     cdir = os.getcwd()
    #     os.chdir('test/test_data/')
    #     command = f'syri -c {coords} -d {delta} -r {ref} -q {qry} -k -F T --prefix test_coords'.split()
    #     p = Popen(command, stdout=PIPE, stderr=PIPE)
    #     out = p.communicate()
    #     print(out)
    #     assert all([f in glob("test_coords*") for f in ["test_coordssyri.out", "test_coordssyri.vcf", "test_coordssyri.summary"]])
    #
    # # END

    def test_bam_input(self):
        from subprocess import Popen, PIPE
        from glob import glob
        import os
        bam = 'out.bam'
        ref = "tair10_chr4.fa.gz"
        qry = "ler_chr4.fa.gz"
        # cdir = os.getcwd()
        # os.chdir('test/test_data/')
        command = f'syri -c {bam} -r {ref} -q {qry}  -F B --prefix test_bam'.split()
        p = Popen(command, stdout=PIPE, stderr=PIPE)
        out = p.communicate()
        print(out)
        assert all([f in glob("test_bam*") for f in ["test_bamsyri.out", "test_bamsyri.vcf", "test_bamsyri.summary"]])
        for f in glob("test_bam*"):
            try:
                os.remove(f)
            except OSError as e:
                if e.errno != 2:    # 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
                    raise
        # cdir = os.getcwd(cdir)
    # END
    def test_paf_input(self):
        from subprocess import Popen, PIPE
        from glob import glob
        paf = 'out.paf'
        ref = "tair10_chr4.fa.gz"
        qry = "ler_chr4.fa.gz"
        prefix = "test_paf"
        # cdir = os.getcwd()
        # os.chdir('test/test_data/')
        command = f'syri -c {paf} -r {ref} -q {qry}  -F P --prefix {prefix}'.split()
        p = Popen(command, stdout=PIPE, stderr=PIPE)
        out = p.communicate()
        print(out)
        assert all([f"{prefix}{f}" in glob(f"{prefix}*") for f in ["syri.out", "syri.vcf", "syri.summary"]])
        for f in glob(f"{prefix}*"):
            try:
                os.remove(f)
            except OSError as e:
                if e.errno != 2:    # 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
                    raise
    # END

if __name__ == '__main__':
    os.chdir('test/test_data/')
    unittest.main(verbosity=3)
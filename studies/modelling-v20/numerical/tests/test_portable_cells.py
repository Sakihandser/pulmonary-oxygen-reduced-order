import importlib.util,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'reproduce.py'
s=importlib.util.spec_from_file_location('rep',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class CellComparison(unittest.TestCase):
    def test_float_vector(self):
        pairs=m.comparable_cells('2.6762431989952593;0.9152717300515847','2.6762431989952593;0.9152717300515846')
        self.assertEqual(len(pairs),2)
        self.assertTrue(all(abs(x-y)<=2e-6+1e-8*abs(x) for x,y in pairs))
    def test_text_changed(self):
        with self.assertRaises(AssertionError):m.comparable_cells('brain;gut','brain;kidney')
    def test_vector_length(self):
        with self.assertRaises(AssertionError):m.comparable_cells('1;2','1;2;3')
    def test_large_numeric_difference(self):
        self.assertFalse(all(abs(x-y)<=2e-6+1e-8*abs(x) for x,y in m.comparable_cells('1;2','1;2.001')))
if __name__=='__main__':unittest.main()

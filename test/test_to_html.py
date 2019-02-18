import unittest
import to_html

class TestToHtml(unittest.TestCase):

    def test_entropy(self):
        self.assertEqual(to_html.entropy([0.25]*4), 2.0)

        self.assertEqual(to_html.entropy([1.0, 0, 0, 0]), 0.0)

        with self.assertRaises(AssertionError):
            to_html.entropy([0.0])

        with self.assertRaises(AssertionError):
            to_html.entropy([2.0, -1.0])

if __name__ == '__main__':
    unittest.main()

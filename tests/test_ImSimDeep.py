"""
Example unit tests for ImSimDeep package
"""
import unittest
import desc.imsimdeep

class ImSimDeepTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'
        
    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.imsimdeep.ImSimDeep(self.message)
        self.assertEquals(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.imsimdeep.ImSimDeep)
        foo = desc.imsimdeep.ImSimDeep(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()

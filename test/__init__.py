from unittest import TestCase

import prsbils

class TestPRSbils(TestCase):
    def test_is_string(self):
        s = funniest.joke()
        self.assertTrue(isinstance(s, basestring))
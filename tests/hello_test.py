import pytest


@pytest.mark.parametrize("input, expected", [("Hello, World!", "Hello, World!")])
def test_hello(input, expected):
    assert input == expected

"""

Collection of unit tests for the AttributePrinter class.

"""

from cobraOps.AttributePrinter import AttributePrinter


class TestAttributePrinter():
    """A collection of tests for the AttributePrinter class.

    """

    def test_str_method(self):
        # Create a new AttributePrinter instance
        instance = AttributePrinter()

        # Add some attributes
        instance.a = "the a attribute"
        instance.b = "the b attribute"

        # Run the __str__ method
        string_representation = instance.__str__()

        # Check that the string contains the attribute information
        assert "a = the a attribute" in string_representation
        assert "b = the b attribute" in string_representation

        # Check that the print method doesn't raise an exception
        print(instance)

    def test_subclass_str_method(self):

        # Define a class that extends the AttributePrinter class
        class myTestClass(AttributePrinter):

            def __init__(self):
                # Add some attributes
                self.a = "the a attribute"
                self.b = "the b attribute"

        # Create a new class instance
        instance = myTestClass()

        # Run the __str__ method
        string_representation = instance.__str__()

        # Check that the string contains the attribute information
        assert "a = the a attribute" in string_representation
        assert "b = the b attribute" in string_representation

        # Check that the print method doesn't raise an exception
        print(instance)


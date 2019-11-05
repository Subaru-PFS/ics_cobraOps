"""

AttributePrinter class.

"""


class AttributePrinter:
    """Class that prints all its instance attributes.

    """

    def __str__(self):
        """Prints the instance attributes.

        Returns
        -------
        str
            A string representation of all the instance attributes.

        """
        attributeList = []

        for attribute in sorted(self.__dict__):
            attributeList.append("%s = %s" % (attribute, getattr(self, attribute)))

        return "\n".join(attributeList)

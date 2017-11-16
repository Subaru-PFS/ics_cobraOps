"""

Attribute printer class.

"""

class AttributePrinter:
    """
    
    Class that prints all its instance attributes.
    
    """
    
    def __str__(self):
        """Prints the instance attributes.
        
        Returns
        -------
        str
            A string representation of the instance attributes.
        
        """
        attributeList = []
        
        for attr in sorted(self.__dict__):
            attributeList.append("{} = {}".format(attr, getattr(self, attr)))
        
        return "\n".join(attributeList)

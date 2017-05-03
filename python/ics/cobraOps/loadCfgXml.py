
import xml.etree.ElementTree as ET


def loadCfgXml(fileName):
    tree = ET.parse(fileName)
    return tree.getroot()

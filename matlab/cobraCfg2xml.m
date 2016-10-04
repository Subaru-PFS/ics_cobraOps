function xml = cobraCfg2xml(CobraConfig, outputfile)

struct4xml.ARM_DATA = CobraConfig.ARM_DATA;

struct2xml(struct4xml,outputfile)

xml = 0;
end
import os
import re
dirList = os.listdir(os.getcwd());

fp = open('blocksTable.h','w');
fp.write('#ifndef GA_BLOCKS_TABLE_H\n#define GA_BLOCKS_TABLE_H\n#include \"blockFunctions.h\"\n\nBlockType BlockTypesTable[] = \n{\n');
fp.close();

fp = open('blockFunctions.h','w');
fp.write('#include \"blocks.h\"\n\n');
fp.close();

p1 = re.compile('(\S+)\.ant');

for f in dirList:
    m = p1.match(f);
    if m:
        print m.group(1)
        os.system('antimonyToFunctions ' + m.group(1) + ' ' + f + ' ' + 'blockFunctions.h blocksTable.h');

fp = open('blocksTable.h','a');
fp.write('\n    {0,0,0,0,0,0,0,0,0,0}\n};\n\n#endif\n');
fp.close();
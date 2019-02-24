import regex
# print regex.sub('amagging{e<=1}',  'found', 'amagingqqqqqqqqqqqqq')
# print dir(regex.match('(amagging){e<=1}', 'amagingqqqqqqqqqqqqq'))
# print regex.match('(amazing){e<=1}', 'amaging').groups()


print regex.sub('(ACCGTGGCCCAGGCGGCCCAG|TAGGCATG){e<=1}', 'found', 'QQQTA+GGCATGQQQQQACCGTGG-CCAGGCGGCCCAGQQQQ')  

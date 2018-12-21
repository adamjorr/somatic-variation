#!/usr/bin/env python3

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The following two lines will be needed in every python script:
from intermine.webservice import Service
import sys

service = Service("https://phytozome.jgi.doe.gov/phytomine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view("name", "briefDescription")

names = ['.'.join(str.split(str.split(str.split(x)[-1], sep = '=')[-1], sep = '.')[0:2]) for x in sys.stdin]

# You can edit the constraint values below
query.add_constraint("organism", "LOOKUP", "E. grandis", code = "A")
query.add_constraint("name", "ONE OF", names, code = "B")

for row in query.rows():
    print(row["name"], row["briefDescription"], sep = '\t')

quit()
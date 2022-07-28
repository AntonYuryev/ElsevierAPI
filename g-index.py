from ElsevierAPI.ScopusAPI.scopus import g_index
from ElsevierAPI import load_api_config

LAST_NAME = 'Yuryev'
FIRST_NAME = 'Anton'

if __name__ == "__main__":
    gi, pub_count, most_cited  = g_index(load_api_config(),LAST_NAME,FIRST_NAME)
    print('%s %s published %d articles' % (FIRST_NAME,LAST_NAME,pub_count))
    print('%s %s\'s g-index is: %d' % (FIRST_NAME,LAST_NAME,gi))
    print('%s %s\'s most cited article:\n%s' % (FIRST_NAME,LAST_NAME,most_cited))


from ..utils import get_auth_token,load_api_config
import urllib.parse, requests
from ..ETM_API.references import DocMine
from scibite_toolkit.scibite_search import SBSRequestBuilder as s


class SBSjson(DocMine):
    @staticmethod
    def dump_fname(search_name): return search_name + '.json'


class SBSapi(s):
    def __init__(self,*args,**kwargs):
        '''
        kwargs:
            api_type - defaults to "search"
        '''
        super().__init__()
        self.APIconfig = dict(args[0]) if args else load_api_config()

        self.set_url(self.APIconfig['SBSurl'])
        self.set_auth_url(self.APIconfig['SBStoken_url'])
        token_address = self.auth_url + "/auth/realms/Scibite/protocol/openid-connect/token"
        token = get_auth_token(token_url=token_address,
                               client_id=self.APIconfig['SBSclientID'], 
                               secret=self.APIconfig['SBSecret'],
                               username=self.APIconfig['SBSuser'],
                               password=self.APIconfig['SBSpassword']
                               )
        
        
        
        self.set_oauth2(self.APIconfig['SBSclientID'], self.APIconfig['SBSecret'])
        print()

        '''
        token_address = self.auth_url + "/auth/realms/Scibite/protocol/openid-connect/token"
        token = get_auth_token(token_url=token_address,
                               user=self.APIconfig['SBSclientID'], 
                               secret=self.APIconfig['SBSecret'])
        req = requests.post(token_address, data={"grant_type": "client_credentials", "client_id": self.APIconfig['SBSclientID'],
                                                 "client_secret": self.APIconfig['SBSecret']},
                            headers={"Content-Type": "application/x-www-form-urlencoded"})
        '''

       
 

    def token_refresh(self):
        s.set_url(self.APIconfig['SBSurl'])
        s.set_auth_url(self.APIconfig['SBStoken_url'])
        s.set_oauth2(self.APIconfig['SBSuser'], self.APIconfig['SBSecret'])
        return s  # This is actually regenerating a new builder request object.

if __name__ == "__main__":
    my_s = SBSapi()
    response_getdocs = s.get_docs(query='abstract ~ INDICATION$D008175 AND DRUG$*',limit=100,markup=True)
    print(response_getdocs['data'])
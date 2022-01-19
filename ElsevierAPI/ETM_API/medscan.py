import subprocess

from textblob.blob import TextBlob

USERDIR = "C:\\Users\\Administrator\\Documents\\MedScan Reader Projects\\"
MEDSCANDIR = 'D:\\MEDSCAN\\DevStandard10.bin\\Standard.bin\\'

PROTEIN = 0
SMALLMOL = 1000000
CELLOBJ = 2000000
COMPLEX = 3000000
CELLPROCESS = 4000000
VIRUS = 5000000
ORGANISM = 7000000
DISEASE = 9000000
CELLTYPE = 10000000
TREATMENT = 13000000
CLINPARAM = 15000000
MEDPROCEDURE = 16000000
FUNCTIONALCLASS = 12000000
UNCLASSIFIED = 31000000
DOMAIN = 45000000

class MedScan:
# uses local MedScan installation to markup text and extrcat relations
    license = str()
    organism = 'Mammal'
    #input_type = 'sentence'
    sentence_maxlen = 3500
    def __init__(self,license,org_name = 'Mammal'):
        self.organism = org_name
        self.license = license
        self.objnames = self.__load_objnames()

    def markup_sentence(self,sentence:str):
        locscan = 'D:\\MEDSCAN\\DevStandard10.bin\\Standard.bin\\locscan'
        sentences = list()
        for i in range (0, len(sentence), self.sentence_maxlen):
            sentences.append(sentence[i:i+self.sentence_maxlen])

        if len(sentences) > 1:
            print('Below sentence is too long and was split for MedScan processing:\n%s' % sentence)

        markup = str()
        for s in sentences:
            completed_process = subprocess.run([locscan, "-t"+self.license, "-W"+USERDIR,"-Q"+self.organism, "sentence:"+s], capture_output=True)
            #print(completed_process.stderr)
            complete_markup = completed_process.stdout.decode('utf-8').strip()
            sentence_start = complete_markup.find('msrc\t', 20)
            sentence_start = 0 if sentence_start < 0 else sentence_start + 5
            concept_scan_markup_pos = complete_markup.find('\r\n<\t#', len(s)+25)
            if concept_scan_markup_pos < 0:
                concept_scan_markup_pos = len(complete_markup)
            markup += complete_markup[sentence_start:concept_scan_markup_pos]+ ' '
        return markup.strip()

    def __load_objnames(self):
        to_return = dict()
        with open(MEDSCANDIR+'xfdata\\'+self.organism+'.ObjectNames.tab', 'r', encoding='utf-8') as f:
            line = f.readline()
            while line:
                dic_row = line.split('\t')
                msid = dic_row[0]
                name = dic_row[1] 
                to_return[msid] = name
                line = f.readline()

        with open(USERDIR+self.organism+'.UserNames.tab', 'r', encoding='utf-8') as f:
            line = f.readline()
            while line:
                dic_row = line.split('\t')
                msid = dic_row[0]
                name = dic_row[1]
                to_return[msid] = name
                line = f.readline()
        return to_return


    def find_concepts(self,paragraph:str) -> dict: 
        blob = TextBlob(paragraph)
        sentences = list(blob.sentences)
        to_return = dict()
        for sentence in sentences:
            medscan_markup = self.markup_sentence(str(sentence))
            range2dict = dict() #{id_range:{id:obj_name}}
            markup_pos = medscan_markup.find('ID{')
            while markup_pos > 0:
                markup_start = markup_pos + 3
                id_end = medscan_markup.find('=',markup_start)
                msids = list(medscan_markup[markup_start:id_end].split(','))
                id_range =  (int(msids[0]) // 1000000) * 1000000
                first_msid = 0 if len(msids) == 1 else 1
                markup_end = medscan_markup.find('}',markup_start+5)
                if markup_end < 0: break #hack for broken markup
                for i in range(first_msid, len(msids)):
                    msid = msids[i]
                    try: obj_name = self.objnames[msid]
                    except KeyError:
                        obj_name = medscan_markup[id_end+1:markup_end]
                        print('"%s" with MedScan ID %s doesn\'t have object name' % (obj_name,msid))     
                    try:                
                        range2dict[id_range][msid] = obj_name
                    except KeyError:
                        range2dict[id_range] = {msid:obj_name}

                markup_pos = medscan_markup.find('ID{',markup_end+1)
            

            to_return[medscan_markup] = range2dict
        
        return to_return # {medscan_markup:{id_range:{id:obj_name}}}, str

    @staticmethod
    def get_drugs(range2dict:dict):
        drugs = []
        for r in [PROTEIN,SMALLMOL,COMPLEX,FUNCTIONALCLASS]: # protein concepts are included for biologics
            try:
                msid_drugs = list(range2dict[r].keys())
                drugs += [x[1] for x in msid_drugs]
            except KeyError: continue
        return drugs if drugs else ['']

    @staticmethod
    def get_diseases(range2dict:dict):
        try:
            msid_diseases = list(range2dict[DISEASE].keys())
            return [x[1] for x in msid_diseases]
        except KeyError: return ['']

    @staticmethod
    def get_concept_type(id_range:int):
        range2type = {PROTEIN:'Protein', SMALLMOL:'Compound',COMPLEX:'Protein',CELLPROCESS:'Cell process',DISEASE:'Disease',FUNCTIONALCLASS:'Protein',
                MEDPROCEDURE:'Medical procedure', CELLTYPE:'Cell type', CLINPARAM:'Clinical parameter', CELLOBJ:'Cell object', ORGANISM:'Organism',
                TREATMENT:'Treatment', UNCLASSIFIED:'Unclassified', DOMAIN:'Domain', VIRUS:'Virus'
                }
        try:
            return range2type[id_range]
        except KeyError:
            print('Name for %d range is not implemented' % id_range)
            return NotImplemented

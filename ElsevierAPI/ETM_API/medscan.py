import subprocess

USERDIR = "C:\\Users\\Administrator\\Documents\\MedScan Reader Projects\\"
MEDSCANDIR = 'D:\\MEDSCAN\\DevStandard10.bin\\Standard.bin\\'

class MedScan:
# uses local MedScan installation to markup text and extrcat relations
    license = str()
    organism = 'Mammal'
    input_type = 'sentence'
    sentence_maxlen = 3500
    def __init__(self,license,org_name = 'Mammal'):
        self.organism = org_name
        self.license = license
        self.objnames = self.__load_objnames()

    def markup(self,sentence:str):
        locscan = 'D:\\MEDSCAN\\DevStandard10.bin\\Standard.bin\\locscan'
        sentences = list()
        for i in range (0, len(sentence), self.sentence_maxlen):
            sentences.append(sentence[i:i+self.sentence_maxlen])

        if len(sentences) > 1:
            print('Below sentence is too long and was split for MedScan processing:\n%s' % sentence)

        markup = ''
        for s in sentences:
            completed_process = subprocess.run([locscan, "-t"+self.license, "-W"+USERDIR,"-Q"+self.organism, self.input_type+":"+s], capture_output=True)
            #print(completed_process.stderr)
            complete_markup = completed_process.stdout.decode('utf-8').strip()
            concept_scan_markup_pos = complete_markup.find('\r\n<\t#')
            if concept_scan_markup_pos > 0:
                markup += complete_markup + ' '
            else:
                markup += complete_markup[:concept_scan_markup_pos]+ ' '
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

    def find_concepts(self,text:str, id_ranges:list):
        medscan_markup = self.markup(text)
        range2dict = dict() #{id_range:{id:obj_name}}
        markup_pos = medscan_markup.find('ID{')
        while markup_pos > 0:
            markup_start = markup_pos + 3
            id_end = medscan_markup.find('=',markup_start)
            msids = list(medscan_markup[markup_start:id_end].split(','))
            markup_end = medscan_markup.find('}',markup_start+5)
            if markup_end < 0: break
            for id_range in id_ranges:
                id_range_max = id_range + 999999
                if int(msids[0]) >= id_range and int(msids[0]) <= id_range_max:
                    first_msid = 0 if len(msids) == 1 else 1
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
        
        return range2dict

    @staticmethod
    def get_drugs(text_concepts:list):
        drugs = []
        for r in [0,1000000,3000000,12000000]:
            try:
                drugs += list(text_concepts[r].values())
            except KeyError: continue
        return drugs if drugs else ['']

    @staticmethod
    def get_diseases(text_concepts:list):
        try:
            return list(text_concepts[9000000].values())
        except KeyError: return ['']


class Regression:
    def __init__(self):
        self.id = 1
        self.name = "BDT"
        self.inputFiles = []
        self.tree = "regNtuplizer/simpleRootTree"
        self.method = "BDT"
        self.tmvaTrainingOptions = []
        self.options = ["!H","!V","BoostType=Grad"]
        self.doErrors = False
        self.doCombine = False
        self.variablesEB = []
        self.variablesEE = []
        self.variablesComb = []
        self.target = "el_genOverSCEnergy"
        self.targetError = "1.4826*abs(BDTresponse - TARGET)"
        self.targetComb = "(el_generatedEnergy-el_scl_rawEnergy*(1+el_scl_preshowerEOverRaw)*BDTresponse)/(el_gsftrk_pAtVtx-el_scl_rawEnergy*(1+el_scl_preshowerEOverRaw)*BDTresponse)"
        self.cuts = []
        self.cutsEB = []
        self.cutsEE = []
        self.cutsError = []
        self.cutsComb = []

    def nameToString(self):
        return "Regression."+str(self.id)+".Name: "+self.name+"\n"

    def inputToString(self):
        if len(self.inputFiles)==0:
            raise StandardError("ERROR: Regression: No input files")
        inputFiles = ""
        for ip in self.inputFiles:
            inputFiles += ip+":"
        inputFiles = inputFiles[0:-1]
        return "Regression."+str(self.id)+".InputFiles: "+inputFiles+"\n"

    def treeToString(self):
        return "Regression."+str(self.id)+".Tree: "+self.tree+"\n"

    def methodToString(self):
        return "Regression."+str(self.id)+".Method: "+self.method+"\n"

    def trainingToString(self):
        if len(self.tmvaTrainingOptions)==0:
            return ""
        options = ""
        for op in self.tmvaTrainingOptions:
            options += op+":"
        options = options[0:-1]
        return "Regression."+str(self.id)+".TMVATrainingOptions: "+options+"\n"

    def optionsToString(self):
        if len(self.options)==0:
            return ""
        options = ""
        for op in self.options:
            options += op+":"
        options = options[0:-1]
        return "Regression."+str(self.id)+".Options: "+options+"\n"

    def variablesEBToString(self):
        if len(self.variablesEB)==0:
            return ""
        var = ""
        for v in self.variablesEB:
            var += v+":"
        var = var[0:-1]
        return "Regression."+str(self.id)+".VariablesEB: "+var+"\n"

    def variablesEEToString(self):
        if len(self.variablesEE)==0:
            return ""
        var = ""
        for v in self.variablesEE:
            var += v+":"
        var = var[0:-1]
        return "Regression."+str(self.id)+".VariablesEE: "+var+"\n"

    def variablesCombToString(self):
        if len(self.variablesComb)==0:
            return ""
        var = ""
        for v in self.variablesComb:
            var += v+":"
        var = var[0:-1]
        return "Regression."+str(self.id)+".VariablesComb: "+var+"\n"

    def targetToString(self):
        return "Regression."+str(self.id)+".Target: "+self.target+"\n"

    def targetErrorToString(self):
        if not self.doErrors:
            return ""
        return "Regression."+str(self.id)+".TargetError: "+self.targetError.replace("TARGET",self.target)+"\n"

    def targetCombToString(self):
        if not self.doCombine:
            return ""
        return "Regression."+str(self.id)+".TargetComb: "+self.targetComb.replace("TARGET",self.target)+"\n"

    def cutsToString(self):
        if len(self.cuts)==0:
            return ""
        cut = ""
        for c in self.cuts:
            cut += c+"&&"
        cut = cut[0:-2]
        return "Regression."+str(self.id)+".CutBase: "+cut+"\n"

    def doErrorsToString(self):
        return "Regression."+str(self.id)+".DoErrors: "+str(self.doErrors)+"\n"

    def doCombineToString(self):
        return "Regression."+str(self.id)+".DoCombine: "+str(self.doCombine)+"\n"

    def cutsEBToString(self):
        if len(self.cutsEB)==0:
            return ""
        cut = ""
        for c in self.cutsEB:
            cut += c+"&&"
        cut = cut[0:-2]
        return "Regression."+str(self.id)+".CutEB: "+cut+"\n"

    def cutsEEToString(self):
        if len(self.cutsEE)==0:
            return ""
        cut = ""
        for c in self.cutsEE:
            cut += c+"&&"
        cut = cut[0:-2]
        return "Regression."+str(self.id)+".CutEE: "+cut+"\n"

    def cutsErrorToString(self):
        if not self.doErrors:
            return ""
        if len(self.cutsError)==0:
            return ""
        cut = ""
        for c in self.cutsError:
            cut += c+"&&"
        cut = cut[0:-2]
        return "Regression."+str(self.id)+".CutError: "+cut+"\n"

    def cutsCombToString(self):
        if not self.doCombine:
            return ""
        if len(self.cutsComb)==0:
            return ""
        cut = ""
        for c in self.cutsComb:
            cut += c+"&&"
        cut = cut[0:-2]
        return "Regression."+str(self.id)+".CutComb: "+cut+"\n"

    def toString(self):
        par = ""
        par += self.nameToString()
        par += self.inputToString()
        par += self.treeToString()
        par += self.methodToString()
        par += self.trainingToString()
        par += self.optionsToString()
        par += self.doErrorsToString()
        par += self.doCombineToString()
        par += self.variablesEBToString()
        par += self.variablesEEToString()
        par += self.variablesCombToString()
        par += self.targetToString()
        par += self.targetErrorToString()
        par += self.targetCombToString()
        par += self.cutsToString()
        par += self.cutsEBToString()
        par += self.cutsEEToString()
        par += self.cutsErrorToString()
        par += self.cutsCombToString()
        return par

    def diff(self, other):
        diff = ""
        if self.id!=other.id:
            diff += "Different ids: "+str(self.id)+" != "+str(other.id)+"\n"
        if self.name!=other.name:
            diff += "Different names: "+self.name+" != "+other.name+"\n"
        if len(self.inputFiles)!=len(other.inputFiles):
            diff += "Different number of input files: "+str(len(self.inputFiles))+" != "+str(len(other.inputFiles))+"\n"
        for file1,file2 in zip(self.inputFiles, other.inputFiles):
            if file1!=file2:
                diff += "Different file:\n"
                diff += "   "+file1+"\n"
                diff += "!= "+file2+"\n"
        if self.tree!=other.tree:
            diff += "Different trees: "+self.tree+" != "+other.tree+"\n"
        diffOptions1 = list(set(self.options) - set(other.options))
        diffOptions2 = list(set(other.options) - set(self.options))
        if len(diffOptions1)>0 or len(diffOptions2)>0:
            diff += "Different options:\n"
            diff += "   "
            for opt in self.options:
                diff += opt+" "
            diff += "\n!= "
            for opt in other.options:
                diff += opt+" "
            diff += "\n"
        diffVarEB1 = list(set(self.variablesEB) - set(other.variablesEB))
        diffVarEB2 = list(set(other.variablesEB) - set(self.variablesEB))
        if len(diffVarEB1)>0 or len(diffVarEB2)>0:
            diff += "Different EB variables:\n"
            if len(diffVarEB1)>0:
                diff += "   New: "
                for var in diffVarEB1:
                    diff += var+" "
                diff += "\n"
            if len(diffVarEB2)>0:
                diff += "   Old: "
                for var in diffVarEB2:
                    diff += var+" "
                diff += "\n"
        diffVarEE1 = list(set(self.variablesEE) - set(other.variablesEE))
        diffVarEE2 = list(set(other.variablesEE) - set(self.variablesEE))
        if len(diffVarEE1)>0 or len(diffVarEE2)>0:
            diff += "Different EE variables:\n"
            if len(diffVarEE1)>0:
                diff += "   New: "
                for var in diffVarEE1:
                    diff += var+" "
                diff += "\n"
            if len(diffVarEE2)>0:
                diff += "   Old: "
                for var in diffVarEE2:
                    diff += var+" "
                diff += "\n"
        diffVarComb1 = list(set(self.variablesComb) - set(other.variablesComb))
        diffVarComb2 = list(set(other.variablesComb) - set(self.variablesComb))
        if len(diffVarComb1)>0 or len(diffVarComb2)>0:
            diff += "Different comb variables:\n"
            if len(diffVarComb1)>0:
                diff += "   New: "
                for var in diffVarComb1:
                    diff += var+" "
                diff += "\n"
            if len(diffVarEE2)>0:
                diff += "   Old: "
                for var in diffVarComb2:
                    diff += var+" "
                diff += "\n"

        if self.target!=other.target:
            diff += "Different targets: "+self.target+" != "+other.target+"\n"
        if self.targetError!=other.targetError:
            diff += "Different error targets: "+self.targetError+" != "+other.targetError+"\n"
        diffCut1 = list(set(self.cuts) - set(other.cuts))
        diffCut2 = list(set(other.cuts) - set(self.cuts))
        if len(diffCut1)>0 or len(diffCut2)>0:
            diff += "Different cuts:\n"
            diff += "   "
            for cut in self.cuts:
                diff += cut+" "
            diff += "\n!= "
            for cut in other.cuts:
                diff += cut+" "
            diff += "\n"
        diffCutVar1 = list(set(self.cutsError) - set(other.cutsError))
        diffCutVar2 = list(set(other.cutsError) - set(self.cutsError))
        if len(diffCut1)>0 or len(diffCut2)>0:
            diff += "Different cutsError:\n"
            diff += "   "
            for cut in self.cutsError:
                diff += cut+" "
            diff += "\n!= "
            for cut in other.cutsError:
                diff += cut+" "
            diff += "\n"
        return diff

import simuOpt, os, sys, time, subprocess, errno,csv
simuOpt.setOptions(alleleType='lineage', quiet=True)
from simuPOP import *
from simuPOP.sampling import drawRandomSample
import datetime
from numpy import *
from itertools import *

homeDir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

options = [
    {'name':'popSize',
     'default':100,
     'label':'Population Size',
     'type': 'integer',
     'validator': 'popSize > 0',
     },
     {'name':'sampSize',
     'default':100,
     'label':'Sample Size',
     'type': 'integer',
     'validator': 'sampSize > 0',
     },
     {'name':'loc',
     'default':5,
     'label':'Number of Loci',
     'type': 'integer',
     'validator': 'loc > 0',
     },
    {'name':'allPerLoc',
     'default':10,
     'label':'Number of Alleles per Loci',
     'type': 'integer',
     'validator': 'allPerLoc > 1',
     },
     {'name':'generations',
     'default':50,
     'label':'Number of Generations',
     'type': 'integer',
     'validator': 'generations >= 0',
     },
     {'name':'replications',
     'default':1000,
     'label':'Number of Replicates',
     'type': 'integer',
     'validator': 'replications > 0',
     },
    {'name': 'dist',
     'default': 'Unif',
     'label': 'Allele Dist.',
     'type': ('chooseOneOf', ['Unif', 'Power', 'Empir']),
     },
    {'name':'runID',
     'default':'runXX',
     'label':'Name of Run',
     'type': 'string',
     'validator': 'isinstance(runID, str)',
     }
   
]

def SaveFstat(pop, output='', maxAllele=0, loci=[], shift=1,
    combine=None):
    if output != '':
        file = output
    else:
        raise exceptions.ValueError, "Please specify output"
    # open file
    try:
        f = open(file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    #
    # file is opened.
    np = pop.numSubPop()
    if np > 200:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 200 samples"
    if loci == []:
        loci = range(pop.totNumLoci())
    nl = len(loci)
    if nl > 100:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 100 loci"
    if maxAllele != 0:
        nu = maxAllele
    else:
        nu = max([max(ind.genotype()) for ind in pop.individuals()]) + 1
    if nu > 999:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 999 alleles at each locus"
        print "If you used simuPOP_la library, you can specify maxAllele in population constructure"
    if nu < 10:
        nd = 1
    elif nu < 100:
        nd = 2
    elif nu < 1000:
        nd = 3
    else: 
        raise ValueError('FSTAT can not handle this now, how many digits?')
    # write the first line
    f.write( '%d %d %d %d\n' % (np, nl, nu, nd) )
    # following lines with loci name.
    for loc in loci:
        f.write( pop.locusName(loc) +"\n");
    gs = pop.totNumLoci()
    for sp in range(0, pop.numSubPop()):
        # genotype of subpopulation sp, individuals are
        # rearranged in perfect order
        gt = pop.genotype(sp)
        for ind in range(0, pop.subPopSize(sp)):
            f.write("%d " % (sp+1))
            p1 = 2*gs*ind          # begining of first hemo copy
            p2 = 2*gs*ind + gs     # second
            for al in loci: # allele
                # change from 0 based allele to 1 based allele
                if combine is None:
                    ale1 = gt[p1+al] + shift
                    ale2 = gt[p2+al] + shift
                    f.write('%%0%dd%%0%dd ' % (nd, nd) % (ale1, ale2))
                else:
                    f.write('%%0%dd' % nd % combine([gt[p1+al], gt[p2+al]]))
            f.write( "\n")
    f.close()
    


def simuGeneticDrift(pathName, popSize=100, sampSize=100, loc=20, allPerLoc = 5, generations=10,  replications=100,  dist = 'unif'):
    alpha =1.5
    if(dist=='unif'):
        chosenFreq = ([1/float(allPerLoc) for number in xrange(allPerLoc)])
    elif(dist=='Power'):
        TotalFreq = sum([pow(number+1,alpha) for number in range(allPerLoc)])
        chosenFreq = [(pow(number+1,alpha)/TotalFreq) for number in range(allPerLoc)]
    else:
        chosenFreq = [0.0113861386, 0.0004950495, 0.01, 0.0381188119, 0.4579207921,  0.2059405941, 0.2049504950, 0.0118811881, 0.0103960396, 0.0500000000 ]
        chosenFreq = chosenFreq/sum(chosenFreq)
        allPerLoc = 10



    lnames = list('Loc' + str(i) for i in range(0,loc))
    pop = Population(size=popSize, loci=loc, ancGen = 0,  lociNames=lnames, infoFields=['ind_id', 'father_id', 'mother_id'])
    pop.dvars().fileName= pathName
    paramList = list([popSize,sampSize,loc,allPerLoc,generations,replications, dist])
    nameList = list(['pop_size','sample_size','n_loci','al_per_loc','n_gen','n_rep', 'allele_dist'])
    paramFile = open(pathName + "/param.csv" ,'wb')
    wr = csv.writer(paramFile, quoting=csv.QUOTE_ALL)
    wr.writerow(nameList)
    wr.writerow(paramList)


    fd = open(pathName + "/runInfo.csv" ,'wb')
    fd.write("rep , file ,gen, r2, Ne1, Ne2, Ne3, Ne4 \n")
    fd.close()

    simu=Simulator(pop, rep=replications)

    pop.setVirtualSplitter(SexSplitter())
    simu.evolve(
        initOps = [
            InitSex(),
            IdTagger(),

            InitGenotype(freq=chosenFreq),

        ],
        preOps = [
            Stat(effectiveSize=ALL_AVAIL, vars='Ne_demo_base_sp',subPops=[0]),
        ],
        postOps = [

            Stat(effectiveSize=ALL_AVAIL, vars='Ne_demo_sp',subPops=[0]),
            Stat(effectiveSize=ALL_AVAIL, LD=list(combinations(range(0, loc), 2)) ,vars=['R2'],subPops=[0]),
            PyOperator(func=preSave, at=-4),
            PyOperator(func=preSave, at=-3),
            PyOperator(func=preSave, at=-2),
            PyOperator(func=coolFun, at=-1)


        ],
        matingScheme = RandomMating(
                ops=[
                    Recombinator(rates = 1),
                    IdTagger(),
                    PedigreeTagger()
                ],
        ),
        gen = generations

   )

    
    return simu


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def coolFun(pop):
    pathName = pop.dvars().fileName
    fd = open(pathName + "/runInfo.csv" ,'a')
    neVal1=(pop.vars()['subPop'][0]['Ne_demo'][0])
    temp = pop.vars()['R2'].items()
    dictlist= list()
    lCount=pop.numLoci()[0]
    for j in range(0,lCount-1):
            dictlist.extend(temp[j][1].values())
    rVal = (mean(dictlist))
    gen = pop.dvars().gen
    rep = pop.dvars().rep

    with open(pathName + "/temp.csv",'r') as temp_file:
        content = [line.rstrip('\n') for line in temp_file]
    temp_file.close()

    myCsvRow = "%s,  sample%s.dat , %s , %s, %s, %s , %s, %s \n" % (rep ,rep ,gen, rVal, neVal1, content[3*rep+2], content[3*rep+1],content[3*rep+0] )
    fd.write(myCsvRow)
    fd.close()




    return(True)

def preSave(pop):
    pathName = pop.dvars().fileName
    fd = open(pathName + "/temp.csv" ,'a')
    neVal1=(pop.vars()['subPop'][0]['Ne_demo'][0])
    myCsvRow = "%s \n" % (neVal1)
    fd.write(myCsvRow)
    fd.close()
    return(True)


if __name__ == '__main__':
    # get all parameters
    pars = simuOpt.Params(options)

    if not pars.getParam():
        sys.exit(0)

    print(pars.asList())
    runName=pars.runID
    pathName = homeDir + '/Samples/' + runName
    print(pathName)
    make_sure_path_exists(pathName)
    sims=simuGeneticDrift(pathName,pars.popSize, pars.sampSize, pars.loc, pars.allPerLoc, pars.generations, pars.replications, pars.dist)



    for i in range(0,pars.replications):
        samp=drawRandomSample(sims.population(i), pars.sampSize)

        outfile = pathName + '/sample' + str(i) + '.dat'

        SaveFstat(samp, outfile)

    os.remove(pathName + "/temp.csv")


    

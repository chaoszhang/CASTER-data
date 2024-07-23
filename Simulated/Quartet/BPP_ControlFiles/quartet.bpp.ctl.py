def generate(file = "bpp.phylip", topID = "ABCD", topology = "(((A,B),C),D);", mbp = 1, mu = 0.01, burnin = 32000, sampfreq = 100, nsample = 4000, thread = 1):
    return """          seed = 2333

       seqfile = """ + file + """
      Imapfile = ../../quartet.Imap.txt
       outfile = """ + file + "." + topID + """.bpp
      mcmcfile = """ + file + "." + topID + """.bpp.mcmc

  speciesdelimitation = 0 * fixed species tree
         speciestree = 1

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  A B C D
                    1 1 1 1
                 """ + topology + """
        phase =   0 0 0 0

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = """ + str(mbp * 1000) + """ * number of data sets in seqfile
         model = gtr

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
* haploid 2uN = 2 * u * 1 -> b = (a - 1) * 2uN
* We make mean = the true value for theta and tau following guides from "A Tutorial on the Use of BPP for Species Tree Estimation and Species Delimitation"
   thetaprior = 3 """ + str(4 * mu) + """ # gamma(a, b) for theta (estimate theta)
     tauprior = 3 """ + str(30.1 * mu) +""" # gamma(a, b) for root tau & Dirichlet(a) for other tau's

      finetune =  1
      
      locusrate = 0 # true to this simulation
      clock = 2 2 2 5 iid g

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = """ + str(burnin) + """
      sampfreq = """ + str(sampfreq) + """
       nsample = """ + str(nsample) + """
       Threads = """ + str(thread)
       
print(generate())

for mu in [0.5, 1, 2]:
    for mbp in [1, 2, 5, 10, 20]:
        with open("data/quartet/" + str(mu) + "X/" + "ABCD" + "_" + str(mbp) + ".bpp.ctl", "w") as f:
            f.write(generate(topID = "ABCD", topology = "(((A,B),C),D);", mbp = mbp, mu = mu * 0.005))
        with open("data/quartet/" + str(mu) + "X/" + "ACBD" + "_" + str(mbp) + ".bpp.ctl", "w") as f:
            f.write(generate(topID = "ACBD", topology = "(((A,C),B),D);", mbp = mbp, mu = mu * 0.005))
        with open("data/quartet/" + str(mu) + "X/" + "BCAD" + "_" + str(mbp) + ".bpp.ctl", "w") as f:
            f.write(generate(topID = "BCAD", topology = "(((B,C),A),D);", mbp = mbp, mu = mu * 0.005))

# locusrate = 1 0 0 5 iid
# the shape parameter (= 5 in the example) specifying how similar u are among loci

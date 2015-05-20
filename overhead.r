#----  Functions to assess a "good" range for the MDL penalty ------------#
#      June 2007, revised Feb. 2009 and Aug. 2010                         #
#                                                                         #
# example of use:  overhead.range(ntax=5) if there are 5 taxa.            #
#                                                                         #
# ntax = number of taxa. No default.                                      #
# nbase = number of bases. Default is 4.                                  #
#     If DNA, this is 4 if gapped sites are excluded and if sites with    #
#     missing data are excluded.                                          #
#     If gapped sites are included, use nbase=5, and run the parsimony    #
#     analysis such that a gap is treated as an additional symbol.        #
#     To include sites missing data: count "?" as an additional symbol    #
#     in the parsimony analysis, and add 1 to nbase.                      #
#-------------------------------------------------------------------------#

log2UB = function(n){
 sum(log2(2*(3:n)-5))
}

penalty.AS2005=function(ntax, nbase=4){
 lgNedge=ceiling(log2(2*ntax-3)); 
 b=ceiling(log2(nbase));
 treedes=2*ntax-4+ntax*ceiling(log2(ntax))
 (treedes+lgNedge)/(b+lgNedge)
}

penalty.integerBits = function(ntax,nbase=4){
 # constrained to integer bits,
 # optimal tree description: by an index < total # trees
 # rather than by the parenthetical description
 lgNedge=ceiling(log2(2*ntax-3)); 
 b=ceiling(log2(nbase));
 treedes = ceiling(log2UB(ntax));
 (treedes+lgNedge)/(b+lgNedge)
}

penalty.optimized = function(ntax,nbase=4){
 # non-integer bits allowed, optimal tree description
 log2UB(ntax+1)/log2(nbase*(2*ntax-3))
}
 
overhead.range = function(ntax,nbase=4){
 cat("Penalty for an additional tree:\n")
 cat(" from AS 2005 paper:              ",penalty.AS2005(ntax,nbase),"\n")
 cat(" with optimized tree description: ",penalty.integerBits(ntax,nbase),"\n")
 cat(" and with non-integer bits:       ",penalty.optimized(ntax,nbase),"\n")
}

import numpy as np

#contributions by Anastasia, Chris, George

#inputs a POSET, outputs a boolean
def hasNumericalSemigroups(P):
    mpvecs = [[a-b for (a,b) in zip(rel[0], rel[1])] for rel in P.MinimalPresentation()]
    if len(mpvecs) == 0:
        return True
    
    T = ToricLattice(len(mpvecs[0]))
    L = T.submodule(mpvecs)
    
    return all(sum([a*b for (a,b) in zip(P.atoms, v)]) % P.m == 0 for v in L.saturation().basis())

# generates those inequalities which Q_{m,b} satisfies
def generateSymmetricPosetInequalities(m,b):
  ineqs = KunzPoset.KunzInequalities(m)
  copy = KunzPoset.KunzInequalities(m)
  for i in inequalities:
    if i[b] == -1:
      copy.append([-1*k for k in i])
  return copy

#returns the cover_relations in the form necessary for the FinitePoset constructor
def maximalSymmetricPoset(m,b):
  Relations = []
  for i in range(m-1):
    if i + 1 != b:
      Relations.append((0, i+1))
      Relations.append((i+1, b))
  return (range(m-1), Relations)

#generates Q_{m,b} as a polyhedron object
def generateSymmetricFace(m,b):
  return Polyhedron(ieqs = generateSymmetricPosetInequalities(m,b))

# inputs m, prints the dimension, f-vector for each Q_{m,b} which b in the range [1, m-1]
def dimensionData(m):
  for b in range(1,m):
    Poly = generateSymmetricFace(m,b)
    print("Q " + str(m) + "," + str(b) + " with dimension " + str(Poly.dimension()) + " and f-vector " + str(Poly.f_vector()))
    
# inputs a face, outputs those inequalities from the Kunz Inequalities which are satisfied with equality.
def facetEqualities(face, m):
  ineqs = KunzPoset.KunzInequalities(m)
  Equalities = []
  for eq in ineqs:
    poly = Polyhedron(eqns = [eq])
    if poly.intersection(face) == facet:
      Equalities.append(eq)
  return Equalities

# Inputs a face and outputs the corresponding KunzPoset - using the KunzPoset constructor
def returnKunzPoset(face,m):
  eqs = FacetEqualities(face,m)
  newEqualities = [iq[1:] for iq in Equalities]
  poset = KunzPoset(m, hyperplane_desc = NewEqualities)
  return poset

# outputs an array of the rays on Q_{m,b} for each b
def findRaysForSymmetryFaces(m):
    KunzCone = Polyhedron(ieqs=KunzPoset.KunzInequalities(m))
    Rays = [tuple(ray) for ray in KunzCone.rays()]
    for b in range(1, m):
        Inequalities = generateSymmetricPosetInequalities(m,b)
        validRays = []
        for ray in Rays:
            check = True
            for ineq in Inequalities:
                if dotProduct(ineq, ray) < 0:
                    check = False
            if check:
                validRays.append(ray)
        print(validRays)
                
# inputs a ray (of C_m), checks all the Kunz Inequalities against the ray, compiles a list of such inequalities and creates a Polyhedron object with those inequalities
  def RayToPolyhedron(ray,m):
    Equalities = []
    for i in range(1,m):
        for j in range(1,m):
            if ray[i-1] + ray[j-1] == ray[(i+j-1) % (m-1)]:
                Eq = [0 for t in range(m)]
                Eq[i]+=1
                Eq[j]+=1
                Eq[((i+j) % (m-1))] =-1
                Equalities.append(Eq)
                Equalities.append([-k for k in Eq])
    return Polyhedron(ieqs = Equalities)
  
  # Input: facet of Q_{m, _} (though can be generalized to any face of C_m) and outputs a dictionary with key: height, value: # of elements in the corresponding Kunz Poset with that height 
  # here height is length of maximal chain from 0 to element
  def posetElementDepth(face, m, showposet = False): # essentially a topological search - definitely not optimal but it doesn't really matter
    KunzP = returnKunzPoset(face.as_polyhedron(), m)
    if showposet:
        KunzP.poset.show()
    Cover_Rel = KunzP.cover_relations
    DepthArray = [-1 for i in range(m)]
    DepthArray[0] = 0
    Depths = {
        0: [0],
        1: [],
        2: [],
        3: [],
        4: []
    }
    CoverRelations = KunzP.cover_relations
    CoverRelDict = {}
    for i in range(0,m):
        CoverRelDict[i] = []
    for rel in Cover_Rel:
        CoverRelDict[rel[0]].append(rel)
        CoverRelDict[rel[1]].append(rel)
    
    currentDepth = 1
    while currentDepth <= 4:
        for x in Depths[currentDepth - 1]:
            for rel in CoverRelDict[x]:
                if rel[0] == x:
                    if rel[1] not in Depths[currentDepth]:
                        Depths[currentDepth].append(rel[1])
                    if DepthArray[rel[1]] != -1 and DepthArray[rel[1]] < currentDepth:
                        Depths[DepthArray[rel[1]]].remove(rel[1]) # x simplifies to currentDepth -1
                    DepthArray[rel[1]] = currentDepth
        currentDepth+=1
    return DepthArray

  # the purpose of the following function is to classify the different types of Kunz Posets obtained from the corresponding facets of Q_{m, _}
  def getTuple(facet, m): 
    DepthArray = PosetElementDepth(facet, m)
    maxDepth = np.argmax(DepthArray) # this is of course the thing at the top, I don't want to count that
    Array = []
    Dict = {
        0: [0],
        1: [],
        2: [],
        3: [],
        4: []
    }
    for x in range(len(DepthArray)):
        #print(DepthArray[x])
        if 1 < DepthArray[x] and DepthArray[x] < DepthArray[maxDepth]:
            Array.append(DepthArray[x])
            #print(str(DepthArray[x]) + " and " + str(x))
            if Dict[DepthArray[x]]:
                Dict[DepthArray[x]].append(x)
            else:
                Dict[DepthArray[x]] = [x]
    ValidDepths = []
    [ValidDepths.append(x) for x in Array if x not in ValidDepths]
    for x in ValidDepths:
        print("Depth " + str(x) + ": " + str(Dict[x]))
    return [Dict[x] for x in ValidDepths]

# applies the previous function for each facet of Q_{m,b}
def getTuples(m,b):
    for facet in generateSymmetryFace(m,b).facets():
        GetTuple(facet, m)

# Easier to read classification of the facet-types of Q_{m,b}. Here a 1-1 boosted element refers to the triple (b/4, b/4, b/2)
def getTypesOfTuples(m,b, printer = False):
    Array = [GetTuple(facet, m) for facet in generateSymmetryFace(m,b).facets()]
    Types = {
        1: [],
        2: [],
        "2B": [],
        3: [],
        "Other": []
    }
    for x in Array:
        if len(x) == 1:
            if len(x[0]) <= 3:
                Types[len(x[0])].append(x)
            else:
                Types["Other"].append(x)
        elif len(x) == 2:
            Types["2B"].append(x)
    if printer:
        print(Types)
        print("Num Facets with 1 boosted element: "+ str(len(Types[1])))
        print("Num Facets with 2 boosted elements: "+ str(len(Types[2])))
        print("Num Facets with 1-1 boosted element: "+ str(len(Types["2B"])))
        print("Num Facets with 3 boosted elements: "+ str(len(Types[3])))
        print("Num Facets with OTHER boosted element: "+ str(len(Types["Other"])))
    return (len(Types[1]), len(Types[2]), len(Types["2B"]), len(Types[3]), len(Types["Other"]))
  
#chris wrote
#inputs a POSET, outputs a boolean
def hasNumericalSemigroups(P):
    mpvecs = [[a-b for (a,b) in zip(rel[0], rel[1])] for rel in P.MinimalPresentation()]
    if len(mpvecs) == 0:
        return True
    
    T = ToricLattice(len(mpvecs[0]))
    L = T.submodule(mpvecs)
    
    return all(sum([a*b for (a,b) in zip(P.atoms, v)]) % P.m == 0 for v in L.saturation().basis())

 #input m, b, return the qmbPolyhedron
#starts with the big list of all inequalities for the kunz cone
#and uses gSPI to make the Qmb inequalities
#pretty sure we can just return qmb
#if shit breaks uncomment the code
#but it just seems redundant
def qmbPolyhedron(m, b):
    bigList=KunzPoset.KunzInequalities(m)
    smallList = generateSymmetricPosetInequalities(m, b)
    qmb = Polyhedron(ieqs = bigList, eqns = smallList)
    return qmb
  
#input m, b, prints posets of the facets of that qmb
#builds the qmb, grabs the facets, and grabs the equalities of the qmb
#now goes thorugh each facet, builds its poset
#also goes through and finds the extra equalities that this facet satisfies
def qmbFacetPosets(m,b):
    poly = qmbPolyhedron(m, b)
    facetList = poly.facets()
    equalitiesList = facetEqualities(m, poly)
    for facet in facetList:
        eqList = facetEqualities(m, facet.as_polyhedron())
        kPoset = KunzPoset(m, hyperplane_desc = eqList)
        #print(facet.as_polyhedron().dimension())
        for eq in eqList:
            if eq not in equalitiesList:
                print(oneEqtoVal(eq))
        kPoset.poset.show()
        
#input m, b, output a list of lists
#the triples allow duplicates, but do not allow the case where 3|b & 3|m and each 3x=3y=3z=b
def findTriples(m, b):
    lst = []
    thirdSet = set()
    halfSet = set()
    if b % 2 == 0:
        halfSet.add((b)/2)
    if m % 2 == b % 2:
            halfSet.add((m+b)/2)
    divByThree = b % 3 == 0 and m % 3 == 0
    if divByThree:
        thirdSet = {b/3, (m+b)/3, (2*m+b)/3}
    for x in [1..m-1]:
        if x !=b: #don't count b
            for y in [x..m-1]:
                if(y != b): #don't count b
                    for z in [y..m-1]:
                        if(z != b and (x+y+z)%m == b): # don't count b and check if they sum to b
                            if  not (x != y and x in thirdSet and y in thirdSet and y in thirdSet): #check if they're the bad kind of triple
                                #now we need to check if one of them is half of b and if the other two are equal
                                if not ((x in halfSet and y != z) or (y in halfSet and x != z) or(z in halfSet and x!=y)):
                                    lst.append([x, y, z])
    return lst
 
#input m, b, output int[][][]
#first finds triples, then for each triple, calls tTE which sends a list of lists with the three equations
#for those three numbers. each of those list[][] get put into one big list and returned
def findTriplesDiffFormat(m,b):
    lst = findTriples(m,b)
    biggest = []
    for triple in lst:
        biggest.append(tripleToEqs(m, triple))
    return biggest
  
#input m, 1-d list, output list[][]
#given the m and a list of three integers, makes the equations in the format that kunzposets and polyhedrons need
#makes it into a tuple so i can have the given length already
#then adds in the x and y and also the sum
def tripleToEqs(m, lst):
    rtn = []
    for x in [0..2]:
        for y in [x+1..2]:
            tup = (0,)*m
            tup = list(tup)
            tup[lst[x]] = 1
            tup[lst[y]] = 1
            tup[((lst[x] + lst[y])%m)] = -1
            rtn.append(tup)
    return rtn
  
#input a 2-d list of equations, outputs a list of the values added
#uses a dictionary to keep track of values, then sees how many times
#we need to add to the list
def eqsToVals(listOfLists):
    vals = {}
    for baby in listOfLists:
        currVals = oneEqtoVal(baby)
        for k,v in currVals.items():
            if k not in vals or v > 1:
                vals[k] = v
    rtn = []
    for k, v in vals.items():
        if v == 2:
            rtn.append(k)
        if len(vals.keys()) == 1:
            rtn.append(k)
        rtn.append(k)
                
                
    return rtn

#input a list of one equation, returns the value(s) that are added in the eq
#inner method of eqToVals
def oneEqtoVal(lst):
    vals = {}
    for i in range(len(lst)):
            if(lst[i] == 1):
                vals[i+1] = 1
            if lst[i] == 2:
                vals[i+1]=2
    return vals
  
#input m, b, output a list of lists, where each inner list is a triple
#representing that facet
#done by building qmb, grabbing facets, then for each facet, finding each
#new equation and finding its values
def facetVals(m, b):
    poly = qmbPolyhedron(m, b)
    facetList = poly.facets()
    qmbEqs = facetEqualities(m, poly)
    rtn = []
    for facet in facetList:
        facetEqualities = facetEqualities(m, facet.as_polyhedron())
        smallerList = []
        for eq in facetEqualities:
            if eq not in qmbEqs:
                smallerList.append(eq)
        vals = eqsToVals(smallerList)
        rtn.append(vals)
    return rtn
  
#input m, output if the facets of all Qmbs with that m have faces
#prints, not returns, uses Chris's semigroups method
def doFacetsHaveSGs(m):
    for b in [1..m-1]:
        poly = qmbPolyhedron(m, b)
        facetList = poly.facets()
        print("The following facets are from Q %d %d" %(m, b))
        for facet in facetList:
            eqList = facetEqualities(m, facet.as_polyhedron())
            kPoset = KunzPoset(m, hyperplane_desc = eqList)
            answer = hasNumericalSemigroups(kPoset)
            print("does this facet have a semigroup?", answer)
            
#input m, b, bPrime, return a list of rays
#makes both polyhedrons, takes the intersection, and grabs the rays of that face
def intersectionRays(m, b, bPrime):
    qmb = qmbPolyhedron(m, b)
    qmbPrime = qmbPolyhedron(m, bPrime)
    qInter = qmb.intersection(qmbPrime)
    return qInter.rays()
  
#input a list of rays, output the list of coordinates that are zero
# throughout all rays, (the subgroup)
def intersectionZeroes(rays):
    listZeroes = []
    if len(rays) == 0:
        return "No intersection"
    for i in range(0, len(rays[0])):
        sum = 0
        for j in range(0, len(rays)):
            sum += rays[j][i]
        if sum == 0:
            listZeroes.append(i+1)
    return listZeroes
  
#input m, b, bPrime, prints the poset of that intersection
#we use intersectionRays and intersectionZeroes to 
def intersectionPoset(m, b, bPrime):
    rays = intersectionRays(m, b, bPrime)
    newrays = []
    firstZero = intersectionZeroes(rays)[0]
    if len(rays) != 0:
        
        for ray in rays:
            newrays.append(ray[0:firstZero-1])
        newM = len(newrays[0]) + 1
        kunzIneqs = KunzPoset.KunzInequalities(newM)
        newEqualities = []
        for i in range(0, len(kunzIneqs)):
            sum = 0
            for j in range(0, len(newrays)):
                sum += dotProduct(kunzIneqs[i],newrays[j])
        #print(sum)
            if sum == 0:
                newEqualities.append(kunzIneqs[i])
        newEqualities=[eq[1:] for eq in newEqualities]
        kPoset = KunzPoset(newM, hyperplane_desc=newEqualities)
        kPoset.poset.show()

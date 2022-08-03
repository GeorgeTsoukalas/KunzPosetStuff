import numpy as np

def generateSymmetricPosetInequalities(m,b): # generates those inequalities which Q_{m,b} satisfies
  ineqs = KunzPoset.KunzInequalities(m)
  copy = KunzPoset.KunzInequalities(m)
  for i in inequalities:
    if i[b] == -1:
      copy.append([-1*k for k in i])
  return copy

def maximalSymmetricPoset(m,b): #returns the cover_relations in the form necessary for the FinitePoset constructor
  Relations = []
  for i in range(m-1):
    if i + 1 != b:
      Relations.append((0, i+1))
      Relations.append((i+1, b))
  return (range(m-1), Relations)

def generateSymmetricFace(m,b): #generates Q_{m,b} as a polyhedron object
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
  def PosetElementDepth(face, m, showposet = False): # essentially a topological search - definitely not optimal but it doesn't really matter
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
  def GetTuple(facet, m): 
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
def GetTuples(m,b):
    for facet in generateSymmetryFace(m,b).facets():
        GetTuple(facet, m)

# Easier to read classification of the facet-types of Q_{m,b}. Here a 1-1 boosted element refers to the triple (b/4, b/4, b/2)
def GetTypesOfTuples(m,b, printer = False):
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

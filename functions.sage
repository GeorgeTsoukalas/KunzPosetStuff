def generateSymmetricPosetInequalities(m,b): # generates those inequalities which Q_{m,b} satisfies
  ineqs = KunzPoset.KunzInequalities(m)
  copy = KunzPoset.KunzInequalities(m)
  for i in inequalities:
    if i[b] == -1:
      copy.append([-1*k for k in i])
  return copy

def MaxSymmetricPoset(m,b): #returns the cover_relations in the form necessary for the FinitePoset constructor
  Relations = []
  for i in range(m-1):
    if i + 1 != b:
      Relations.append((0, i+1))
      Relations.append((i+1, b))
  return (range(m-1), Relations)

def generateSymmetricFace(m,b): #generates Q_{m,b} as a polyhedron object
  return Polyhedron(ieqs = generateSymmetricPosetInequalities(m,b))

def dimensionData(m):
  for b in range(1,m):
    Poly = generateSymmetricFace(m,b)
    print("Q " + str(m) + "," + str(b) + " with dimension " + str(Poly.dimension()) + " and f-vector " + str(Poly.f_vector()))
    
def facetEqualities(face, m):
  ineqs = KunzPoset.KunzInequalities(m)
  Equalities = []
  for eq in ineqs:
    poly = Polyhedron(eqns = [eq])
    if poly.intersection(face) == facet:
      Equalities.append(eq)
  return Equalities

def returnKunzPoset(face,m):
  eqs = FacetEqualities(face,m)
  newEqualities = [iq[1:] for iq in Equalities]
  poset = KunzPoset(m, hyperplane_desc = NewEqualities)
  return poset

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
  
  def PosetElementDepth(facet, m, showposet = False): # essentially a topological search
    KunzP = returnKunzPoset(facet.as_polyhedron(), m)
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

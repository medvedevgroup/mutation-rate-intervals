#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from gzip import open as gzip_open


def usage(s=None):
	message = """
usage: cat <simulation_table> | evaluate_hypergeometric_slicer [options]
  --confidence=<p>          (C=) size of confidence interval; if present, this
                            overrides the value in the input file
                            (default is to get this from the input file)
  --slices=<N>              (m=) number of slices
                            (default is 100)
  --maxsketch=<N>           maximum sketch size; note that sketch sizes are
                            defined in the input file; this option causes us
                            to ignore large sketches
  --useL.A,B                use the column named "L.A,B" instead of L; input is
                            still required to contain an "L" column; the L.A,B
                            value is used as L in computing the confidence
                            interval but the original L is reported as L in the
                            output; L.A,B is reported as L in each record of
                            the details file
  --details=<filename>      write record-by-record details to a file
                            (by default we do not report these)
  --progress=<number>       periodically report how many input records we've
                            processed

typical input:
  #L      K  r     confidence q           nIntersection(s=100) nIntersection(s=500) ...
  4500000 21 0.100 0.95       0.890581011 7                    31                   ...
  4500000 21 0.100 0.95       0.890581011 4                    29                   ...
  4500000 21 0.100 0.95       0.890581011 7                    35                   ...
  4500000 21 0.100 0.95       0.890581011 5                    27                   ...
   ...

Columns L, k, r1, confidence, and q are required. At least one nIntersection
column is required. Sketch sizes are inferred fron the nIntersection column
headers.

nIntersection is the number of kmers (or kspans) that are in the intersection
of BS(A), BS(B), and BS(A union B) for a sketch of the given size. The Jaccard
estimate is nIntersection/s."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	confidence         = None
	numSlices          = 100
	useAlternateL      = False
	maxSketchSize      = None
	desiredSketchSizes = None
	whichSlicer        = None
	detailsFilename    = None
	reportProgress     = None
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--confidence=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--slices=")) \
		 or (arg.startswith("m=")) or (arg.startswith("--m=")) \
		 or (arg.startswith("M=")) or (arg.startswith("--M=")):
			numSlices = int(argVal)
		elif (arg.lower() == "--usel.a,b"):
			useAlternateL = True
		elif (arg.startswith("--maxsketch=")):
			maxSketchSize = int_with_unit(argVal)
		elif (arg == "--slicer=standard"):
			whichSlicer = "standard"
		elif (arg == "--slicer=zetamatic"):
			whichSlicer = "zetamatic"
		elif (arg.startswith("--details=")):
			detailsFilename = argVal
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# decide which slicer module to use

	if (whichSlicer in ["standard",None]):
		import hypergeometric_slicer as hgslicer
	elif (whichSlicer == "zetamatic"):
		import hypergeometric_slicer_zetamatic as hgslicer
	else:
		assert (False), "unkown slicer: %s" % whichSlicer

	if ("nocache" in debug):
		hgslicer.useCache = False

	if ("nosanity" in debug):
		hgslicer.useNLowSanityCheck = False
		hgslicer.useNHighSanityCheck = False

	if ("showzeta" in debug):
		hgslicer.showZetaCalls = True

	slicerName = hgslicer.moduleName
	if ("noshortcut" in debug): slicerName += ",noshortcut"

	# process the simulation table

	detailsF = None
	if (detailsFilename != None):
		if (detailsFilename.endswith(".gz")) or (detailsFilename.endswith(".gzip")):
			detailsF = gzip_open(detailsFilename,"wt")
		else:
			detailsF = open(detailsFilename,"wt")

	if (detailsF != None):
		print("#%s" % "\t".join(["trial","r1","k","L","confidence","s","m","q","jHat","r1Left","r1Right","inCI"]),
		      file=detailsF)

	paramsToTrials    = {}
	paramsToSuccesses = {}

	recordNum = 0
	for trial in read_simulation_records(stdin,
	                                     useAlternateL=useAlternateL,
	                                     confidenceOverride=confidence):
		recordNum += 1
		if (reportProgress != None):
			if (recordNum <= 2) or (recordNum % reportProgress == 0):
				print("processing record %d" % recordNum,file=stderr)

		(L,k,r1,confidence,q) = (trial.L,trial.k,trial.r1,trial.confidence,trial.q)
		alpha = 1 - confidence

		if (useAlternateL):
			alternateL = trial.alternateL
			LforCI = alternateL
			strL = "%.1f" % alternateL
		else:
			LforCI = L
			strL = "%d" % L

		for s in trial.nIntersection:
			if (maxSketchSize != None) and (s > maxSketchSize):
				if (reportProgress != None) and (recordNum == 1):
					print("processing record %d, ignoring sketch size %d" % (recordNum,s),file=stderr)
				continue

			if (reportProgress != None) and (recordNum == 1):
				print("processing record %d, sketch size %d" % (recordNum,s),file=stderr)
			params = (L,k,r1,confidence,q,s)
			if (params not in paramsToTrials):
				paramsToTrials[params] = paramsToSuccesses[params] = 0
			jaccardObserved = trial.nIntersection[s] / s

			r1Left = r1Right = float("nan")
			paramsToTrials[params] += 1
			if ("noshortcut" in debug):
				(r1Left,r1Right) = hgslicer.r1_confidence_interval(LforCI,k,s,alpha,numSlices,jaccardObserved)
				success = 1 if (r1Left <= r1 <= r1Right) else 0
				if ("showcalls" in debug):
					print("hgslicer.r1_confidence_interval(%s,%s,%s,%s,%s,%s) = %d" \
					    % (LforCI,k,s,alpha,numSlices,jaccardObserved,success),
					       file=stderr)
			else:
				success = hgslicer.truth_in_jaccard_bounds(LforCI,k,r1,s,alpha,numSlices,jaccardObserved)
				if ("showcalls" in debug):
					print("hgslicer.truth_in_jaccard_bounds(%s,%s,%s,%s,%s,%s,%s) = %d" \
					    % (LforCI,k,r1,s,alpha,numSlices,jaccardObserved,success),
					       file=stderr)
			paramsToSuccesses[params] += success

			if (detailsF != None):
				print("%d %.3f %d %s %.3f %d %d %.9f %.9f %.9f %.9f %d" \
					% (recordNum,r1,k,strL,confidence,s,numSlices,q,jaccardObserved,r1Left,r1Right,success),
					  file=detailsF)

	parameterSets = list(paramsToTrials.keys())
	parameterSets.sort()

	if (detailsF != None):
		detailsF.close()

	print("#%s" % "\t".join(["module","r1","k","L","confidence","s","m","trials","q","inCI"]))
	for params in parameterSets:
		(L,k,r1,confidence,q,s) = params
		print("%s\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.9f\t%.3f" \
		    % (slicerName,
		       r1,k,L,confidence,s,numSlices,
		       paramsToTrials[params],q,paramsToSuccesses[params]/paramsToTrials[params]))


# read_simulation_records--
#	Yield input records one at a time. Typical input is as shown in this
#	program's usage statement.

class Trial: pass

def read_simulation_records(f,useAlternateL=False,confidenceOverride=None):
	headerLine    = None
	numColumnsNeeded = None
	columnNames   = None
	sketchSizes   = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("#")):
			if (headerLine == None):
				headerLine = line
				headerFields = line.split()
				headerLineNumber = lineNumber
				if (columnNames == None):
					fields = line.split()
					fields[0] = fields[0][1:]
					extraRequired = None if (not useAlternateL)        else ["L.A,B"]
					notRequired   = None if (confidenceOverride==None) else ["confidence"]
					(columnNames,sketchSizes) \
					    = decipher_column_names(fields,extraRequired=extraRequired,notRequired=notRequired)
			else:
				assert (line.split() == headerFields), \
				       "inconsistent headers at lines %d and %d" \
				     % (headerLineNumber,lineNumber)
			continue

		assert (columnNames != None), \
		       "input column names are not provided within the input file"

		if (numColumnsNeeded == None):
			numColumnsNeeded = 1 + max([columnNames[name] for name in columnNames])

		fields = line.split()
		assert (len(fields) >= numColumnsNeeded), \
		       "not enough columns at line %d (%d, expected %d)" \
		     % (lineNumber,len(fields),numColumnsNeeded)

		trial = Trial()
		trial.L  = int  (fields[columnNames["L"]])
		trial.k  = int  (fields[columnNames["k"]])
		trial.r1 = float(fields[columnNames["r1"]])
		trial.q  = float(fields[columnNames["q"]])

		if (confidenceOverride == None):
			trial.confidence = float(fields[columnNames["confidence"]])
		else:
			trial.confidence = confidenceOverride

		if (useAlternateL): trial.alternateL = float(fields[columnNames["L.A,B"]])

		trial.nIntersection = {}
		for s in sketchSizes:
			trial.nIntersection[s] = int(fields[columnNames["nIntersection(s=%d)"%s]])

		yield trial


# decipher_column_names--
#	set up a hash from variable name to (zero-based) column number

requiredColumns = ["L","k","r1","confidence","q"]
columnAliases   = {"K" : "k",
                   "r" : "r1"}

def decipher_column_names(names,extraRequired=None,notRequired=None):
	if (extraRequired == None): extraRequired = []
	if (notRequired   == None): notRequired   = []
	columnNames = {}
	sketchSizes = []
	for (ix,name) in enumerate(names):
		actualName = name
		if (name in columnAliases): name = columnAliases[name]
		if (name not in requiredColumns+extraRequired): continue
		if (name in columnNames):
			exit("column name \"%s\" (or an alias) appears more than once" % name)
		columnNames[name] = ix
	for name in requiredColumns+extraRequired:
		if (name in notRequired): continue
		if (name not in columnNames):
			exit("input file lacks required name \"%s\"" % name)
	for (ix,name) in enumerate(names):
		if (name in requiredColumns+extraRequired): continue
		if (name in columnNames):
			exit("column name \"%s\" appears more than once" % name)
		prefix = "nIntersection(s="
		suffix = ")"
		if (name.startswith(prefix)) and (name.endswith(suffix)):
			sketchSizes += [int(name[len(prefix):-len(suffix)])]
		columnNames[name] = ix
	if (sketchSizes == []):
		exit("input file lacks any sketch results")
	return (columnNames,sketchSizes)


# parse_probability--
#	Parse a string as a probability

def parse_probability(s,strict=True):
	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try:
		p = float(s)
	except:
		try:
			(numer,denom) = s.split("/",1)
			p = float(numer)/float(denom)
		except:
			raise ValueError

	p *= scale

	if (strict) and (not 0.0 <= p <= 1.0):
		raise ValueError

	return p


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))

if __name__ == "__main__": main()

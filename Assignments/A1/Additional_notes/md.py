import sys

if __name__ == '__main__': 
	in_filename  = sys.argv[1] if len(sys.argv)>1 else "Some Notes on Assignment One.txt"
	theme        = "Simplex"
	md_string    = open(in_filename).read()
	out_filename = sys.argv[2] if len(sys.argv)>2 else "Some Notes on Assignment One.html"
	open(out_filename,"w").write(\
	"""<!DOCTYPE html>
	<html>
	<title>Some Notes on Assignment One</title>
	<xmp id='md' theme='"""+theme+"""'>"""+md_string+"""</xmp>
	<script src='http://strapdownjs.com/v/0.2/strapdown.js'></script>
	</html>""")

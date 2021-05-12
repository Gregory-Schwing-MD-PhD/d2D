<!DOCTYPE HTML PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xml:lang="EN" xmlns="http://www.w3.org/1999/xhtml" lang="EN"><head>

<title>The Vendruscolo Laboratory</title><meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">

<link rel="stylesheet" href="../vendruscolo.css" type="text/css">
<style type="text/css">
<script language="JavaScript" type="text/javascript">

</style><!-- base href="http://www-vendruscolo.ch.cam.ac.uk/index.html" --></head><body>

<div id="outer">
<div id="inner">

<!--TITLE-->

<div id="banner" align="left">
		  
<a href="http://www.cam.ac.uk/">
<img src="../Images/cambridge_logo2.png" alt="">
</a>
<h1>Carlo Camilloni - The Vendruscolo Laboratory</h1>
<h2> 
<a href="http://www.ch.cam.ac.uk/" style="text-decoration:none">Department of Chemistry </a>
<br>
<a href="http://www.cam.ac.uk/" style="text-decoration:none">University of Cambridge </a>
</h2>
</div>

<!--TOP NAVIGATION BAR-->
<div id="header">
		 <ul id="topMenu">
		 		 <li class="first"><a href="">&nbsp;</a></li>
				 <li><a href="http://www-vendruscolo.ch.cam.ac.uk/index.html">Home</a></li>
				 <li><a href="http://www-vendruscolo.ch.cam.ac.uk/research.html">Research</a></li>
				 <li><a href="http://www-vendruscolo.ch.cam.ac.uk/publications.html">Publications</a></li>
				 <li><a href="http://www-vendruscolo.ch.cam.ac.uk/people.html">People</a></li>
				 <li><a href="http://www-vendruscolo.ch.cam.ac.uk/software.html">Software</a></li>
				 <li><a href="http://www-vendruscolo.ch.cam.ac.uk/contact.html">Contacts</a></li>
				 <li class="last"><a href="">&nbsp;</a></li>	 
		 </ul>
</div>

<div id="body" align="left">
<br>

<p><h2>&delta;2D - Determination of secondary structure populations from chemical shifts</h2></p>
<br>
<p> INPUT file: the input file is a text file following the SHIFTY format </p>

<form enctype="multipart/form-data" action="d2D.cc536.php" method="POST">
<!-- MAX_FILE_SIZE must precede the file input field -->
<input type="hidden" name="MAX_FILE_SIZE" value="7000000" />
<!-- Name of input element determines name in $_FILES array -->
<input name="userfile" type="file"/>
<input type=checkbox name="ph" value="Yes"/>low pH<br><br>
<input type="submit" value="Send File"/>
</form>

<br></br>
<p>EXAMPLE:</p>
<pre>#NUM	AA	HA	CA	CB	CO	N	HN   
1	P	4.32	62.69	32.89	0	0	0    
2	N	5.03	52.39	38.99	174.40	0	0    
3	F	4.32	60.09	31.99	176.00	121.78	9.82 
4	S	4.30	60.39	64.19	174.30	112.38	8.49 
5	G	3.96	44.59	0	170.40	110.48	9.22 
6	N	5.52	52.69	39.89	174.00	118.48	8.15 
7	E	5.14	56.39	31.89	175.00	123.18	9.35 
8	K	5.28	53.79	35.79	174.90	123.08	10.17
9	I	3.92	55.19	33.29	175.80	125.48	9.05 </pre>
<br>
<p> NOTES: </p> 
<p>*Non conventional residue notation: X=Oxidized cysteines.</p>
<p>**Possible pHs: acid (<=3) or neutral (Default: neutral).</p><br/>

<table>
<tr>
<td>
Tables
</td>
</tr>
<tr>
<td>
<p><a href="docs/average-tables.pdf">Average chemical shifts</a></p>
</td>
</tr>
<tr>
<td>
<p><a href="docs/CA-helix-corrections.pdf">CA helix corrections</a></p>
<p><a href="docs/CB-helix-corrections.pdf">CB helix corrections</a></p>
<p><a href="docs/CO-helix-corrections.pdf">CO helix corrections</a></p>
<p><a href="docs/HA-helix-corrections.pdf">HA helix corrections</a></p>
<p><a href="docs/HN-helix-corrections.pdf">HN helix corrections</a></p>
<p><a href="docs/NH-helix-corrections.pdf">N helix corrections</a></p>
</td>
<td>
<p><a href="docs/CA-beta-corrections.pdf">CA beta corrections</a></p>
<p><a href="docs/CB-beta-corrections.pdf">CB beta corrections</a></p>
<p><a href="docs/CO-beta-corrections.pdf">CO beta corrections</a></p>
<p><a href="docs/HA-beta-corrections.pdf">HA beta corrections</a></p>
<p><a href="docs/HN-beta-corrections.pdf">HN beta corrections</a></p>
<p><a href="docs/NH-beta-corrections.pdf">N beta corrections</a></p>
</td>
<td>
<p><a href="docs/CA-coil-corrections.pdf">CA coil corrections</a></p>
<p><a href="docs/CB-coil-corrections.pdf">CB coil corrections</a></p>
<p><a href="docs/CO-coil-corrections.pdf">CO coil corrections</a></p>
<p><a href="docs/HA-coil-corrections.pdf">HA coil corrections</a></p>
<p><a href="docs/HN-coil-corrections.pdf">HN coil corrections</a></p>
<p><a href="docs/NH-coil-corrections.pdf">N coil corrections</a></p>
</td>
<td>
<p><a href="docs/CA-ppii-corrections.pdf">CA PPII corrections</a></p>
<p><a href="docs/CB-ppii-corrections.pdf">CB PPII corrections</a></p>
<p><a href="docs/CO-ppii-corrections.pdf">CO PPII corrections</a></p>
<p><a href="docs/HA-ppii-corrections.pdf">HA PPII corrections</a></p>
<p><a href="docs/HN-ppii-corrections.pdf">HN PPII corrections</a></p>
<p><a href="docs/NH-ppii-corrections.pdf">N PPII corrections</a></p>
</td>
</tr>
</table>
<br>
<p><h3>Reference</h3></p>
<p>C. Camilloni, A. De Simone, W. Vranken and M. Vendruscolo.</p>
<br></br>
<!--<p>Unique accesses since 24 Sep 2010</p>
<script language="Javascript" src="http://www-vendruscolo.ch.cam.ac.uk/d2D/counter.php?page=index"></script>
-->
</div>
</body>
</html>

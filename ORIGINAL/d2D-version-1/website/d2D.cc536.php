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

<?php
	/* In PHP versions earlier than 4.1.0, $HTTP_POST_FILES should be used instead of $_FILES */
	echo '<pre>';
	$uploaddir = 'uploads/';
	$uploadfile = $uploaddir . basename($_FILES['userfile']['name']) . '.tmpx';
	$num=mt_rand();

        if(isset($_POST['ph']) && $_POST['ph'] == 'Yes') 
        {
           echo "Low pH selected";
           $command="./program/d2D.x -file uploads/tmp.$num -pH acid -ppii -out uploads/SS-$num.dat > $uploadfile" . ".out 2> uploads/err.$num";
        } else {
           $command="./program/d2D.x -file uploads/tmp.$num -pH neutral -ppii -out uploads/SS-$num.dat > $uploadfile" . ".out 2> uploads/err.$num";
        }

	echo '<pre>';

	if (move_uploaded_file($_FILES['userfile']['tmp_name'], $uploadfile))
	{
		$l=system("./sy2dd.sh $uploadfile $num");
		$l=system($command);
	} 
	else
	{
		echo "File failed to upload correctly!\n";
	}

	if (@readfile("uploads/SS-$num.dat"))
        { 
          echo "Done!\n\nNOTES:\n";
	  $l=system("grep WARNING uploads/err.$num > uploads/tmp");
	  @readfile("uploads/tmp");
        }
        else
        {
          echo "Error, Check your input file!\n";
	  @readfile("uploads/err.$num");
	}
	$l=system("rm -f uploads/*.$num $uploadfile");
	echo '</pre>';
	echo "<p><hr>Comments and bugs to <a href=\"mailto:cc536@cam.ac.uk\">Carlo Camilloni</a>.<hr></p>"; 
?>

<h2><a href="http://www-vendruscolo.ch.cam.ac.uk/d2D/index.php">Back to d2D</a></h2>
</div>
</body>
</html>

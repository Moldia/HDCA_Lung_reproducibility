setBatchMode(true);

/*set the number horizontal (x) tiles & vertical tiles (y)(IMPROTANT Note: these numbers should be at least "2" and when they divide the corresponding dimensions, the result has to be INTEGER)*/
x=2 
y=2

dir=getDirectory("Choose a Directory"); 
print(dir); 
list = getFileList(dir); 
for (i=0; i<list.length; i++) {
	open(dir+list[i]); 
     title = getTitle;
     print(title);
     title = replace(title, ".tif", "");
	splitDir= dir + title + "\\";
	print(splitDir); 
	File.makeDirectory(splitDir);
	if (endsWith(list[i], ".tif")){ 
               print(i + ": " + dir+list[i]);  
               imgName=getTitle(); 
         baseNameEnd=indexOf(imgName, ".tif"); 
         baseName=substring(imgName, 0, baseNameEnd);
         run("Montage to Stack...", "columns=x rows=y border=0");
			run("Stack to Images");
			selectImage(imgName);
			close();
for (j=0;j<x*y;j++) {
	selectImage(j+1);
     title1 = getTitle;
     print(title1);
     saveAs("tif", splitDir+title1);
	}
  }
  run("Close All");
}


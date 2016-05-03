Clazz.declarePackage ("J.io");
Clazz.load (["J.api.JmolZipUtilities"], "J.io.JmolUtil", ["java.io.BufferedInputStream", "$.BufferedReader", "java.net.URL", "java.util.Hashtable", "$.StringTokenizer", "JU.AU", "$.Lst", "$.OC", "$.PT", "$.Rdr", "$.SB", "J.adapter.smarter.AtomSetCollection", "J.api.Interface", "J.io.JmolBinary", "JU.Escape", "$.Logger"], function () {
c$ = Clazz.declareType (J.io, "JmolUtil", null, J.api.JmolZipUtilities);
Clazz.makeConstructor (c$, 
function () {
});
c$.checkSpecialData = Clazz.defineMethod (c$, "checkSpecialData", 
 function (zpt, is, zipDirectory) {
var isSpartan = false;
for (var i = 1; i < zipDirectory.length; i++) {
if (zipDirectory[i].endsWith (".spardir/") || zipDirectory[i].indexOf ("_spartandir") >= 0) {
isSpartan = true;
break;
}}
if (!isSpartan) return null;
var data =  new JU.SB ();
data.append ("Zip File Directory: ").append ("\n").append (JU.Escape.eAS (zipDirectory, true)).append ("\n");
var fileData =  new java.util.Hashtable ();
zpt.getAllZipData (is,  Clazz.newArray (-1, []), "", "Molecule", fileData);
var prefix = "|";
var outputData = fileData.get (prefix + "output");
if (outputData == null) outputData = fileData.get ((prefix = "|" + zipDirectory[1]) + "output");
data.append (outputData);
var files = J.io.JmolUtil.getSpartanFileList (prefix, J.io.JmolUtil.getSpartanDirs (outputData));
for (var i = 2; i < files.length; i++) {
var name = files[i];
if (fileData.containsKey (name)) data.append (fileData.get (name));
 else data.append (name + "\n");
}
return data;
}, "javajs.api.GenericZipTools,java.io.InputStream,~A");
c$.checkSpecialInZip = Clazz.defineMethod (c$, "checkSpecialInZip", 
function (zipDirectory) {
var name;
return (zipDirectory.length < 2 ? null : (name = zipDirectory[1]).endsWith (".spardir/") || zipDirectory.length == 2 ?  Clazz.newArray (-1, ["", (name.endsWith ("/") ? name.substring (0, name.length - 1) : name)]) : null);
}, "~A");
c$.getSpartanDirs = Clazz.defineMethod (c$, "getSpartanDirs", 
 function (outputFileData) {
if (outputFileData == null) return  Clazz.newArray (-1, []);
var v =  new JU.Lst ();
var token;
var lasttoken = "";
if (!outputFileData.startsWith ("java.io.FileNotFoundException") && !outputFileData.startsWith ("FILE NOT FOUND") && outputFileData.indexOf ("<html") < 0) try {
var tokens =  new java.util.StringTokenizer (outputFileData, " \t\r\n");
while (tokens.hasMoreTokens ()) {
if ((token = tokens.nextToken ()).equals (")")) v.addLast (lasttoken);
 else if (token.equals ("Start-") && tokens.nextToken ().equals ("Molecule")) v.addLast (JU.PT.split (tokens.nextToken (), "\"")[1]);
lasttoken = token;
}
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
} else {
throw e;
}
}
return (v.size () == 0 ?  Clazz.newArray (-1, ["M0001"]) : v.toArray ( new Array (v.size ())));
}, "~S");
c$.getSpartanFileList = Clazz.defineMethod (c$, "getSpartanFileList", 
 function (name, dirNums) {
var files =  new Array (2 + dirNums.length * 5);
files[0] = "SpartanSmol";
files[1] = "Directory Entry ";
var pt = 2;
name = name.$replace ('\\', '/');
if (name.endsWith ("/")) name = name.substring (0, name.length - 1);
var sep = (name.endsWith (".zip") ? "|" : "/");
for (var i = 0; i < dirNums.length; i++) {
var path = name + sep;
path += (JU.PT.isDigit (dirNums[i].charAt (0)) ? "Profile." + dirNums[i] : dirNums[i]) + "/";
files[pt++] = path + "#JMOL_MODEL " + dirNums[i];
files[pt++] = path + "input";
files[pt++] = path + "archive";
files[pt++] = path + "Molecule:asBinaryString";
files[pt++] = path + "proparc";
}
return files;
}, "~S,~A");
c$.shortSceneFilename = Clazz.defineMethod (c$, "shortSceneFilename", 
 function (pathName) {
var pt = pathName.indexOf ("_scene_") + 7;
if (pt < 7) return pathName;
var s = "";
if (pathName.endsWith ("|state.spt")) {
var pt1 = pathName.indexOf ('.', pt);
if (pt1 < 0) return pathName;
s = pathName.substring (pt, pt1);
}var pt2 = pathName.lastIndexOf ("|");
return pathName.substring (0, pt) + s + (pt2 > 0 ? pathName.substring (pt2) : "");
}, "~S");
Clazz.overrideMethod (c$, "getAtomSetCollectionOrBufferedReaderFromZip", 
function (vwr, adapter, is, fileName, zipDirectory, htParams, subFilePtr, asBufferedReader) {
var doCombine = (subFilePtr == 1);
htParams.put ("zipSet", fileName);
var subFileList = htParams.get ("subFileList");
if (subFileList == null) subFileList = J.io.JmolUtil.checkSpecialInZip (zipDirectory);
var subFileName = (subFileList == null || subFilePtr >= subFileList.length ? null : subFileList[subFilePtr]);
if (subFileName != null && (subFileName.startsWith ("/") || subFileName.startsWith ("\\"))) subFileName = subFileName.substring (1);
var selectedFile = 0;
if (subFileName == null && htParams.containsKey ("modelNumber")) {
selectedFile = (htParams.get ("modelNumber")).intValue ();
if (selectedFile > 0 && doCombine) htParams.remove ("modelNumber");
}var manifest = htParams.get ("manifest");
var useFileManifest = (manifest == null);
if (useFileManifest) manifest = (zipDirectory.length > 0 ? zipDirectory[0] : "");
var haveManifest = (manifest.length > 0);
if (haveManifest) {
if (JU.Logger.debugging) JU.Logger.debug ("manifest for  " + fileName + ":\n" + manifest);
}var ignoreErrors = (manifest.indexOf ("IGNORE_ERRORS") >= 0);
var selectAll = (manifest.indexOf ("IGNORE_MANIFEST") >= 0);
var exceptFiles = (manifest.indexOf ("EXCEPT_FILES") >= 0);
if (selectAll || subFileName != null) haveManifest = false;
if (useFileManifest && haveManifest) {
var path = J.io.JmolBinary.getManifestScriptPath (manifest);
if (path != null) return "NOTE: file recognized as a script file: " + fileName + path + "\n";
}var vCollections =  new JU.Lst ();
var htCollections = (haveManifest ?  new java.util.Hashtable () : null);
var nFiles = 0;
var zpt = vwr.getJzt ();
var ret = J.io.JmolUtil.checkSpecialData (zpt, is, zipDirectory);
if (Clazz.instanceOf (ret, String)) return ret;
var data = ret;
try {
if (data != null) {
var reader = JU.Rdr.getBR (data.toString ());
if (asBufferedReader) return reader;
ret = adapter.getAtomSetCollectionFromReader (fileName, reader, htParams);
if (Clazz.instanceOf (ret, String)) return ret;
if (Clazz.instanceOf (ret, J.adapter.smarter.AtomSetCollection)) {
var atomSetCollection = ret;
if (atomSetCollection.errorMessage != null) {
if (ignoreErrors) return null;
return atomSetCollection.errorMessage;
}return atomSetCollection;
}if (ignoreErrors) return null;
return "unknown reader error";
}if (Clazz.instanceOf (is, java.io.BufferedInputStream)) is = JU.Rdr.getPngZipStream (is, true);
var zis = zpt.newZipInputStream (is);
var ze;
if (haveManifest) manifest = '|' + manifest.$replace ('\r', '|').$replace ('\n', '|') + '|';
while ((ze = zis.getNextEntry ()) != null && (selectedFile <= 0 || vCollections.size () < selectedFile)) {
if (ze.isDirectory ()) continue;
var thisEntry = ze.getName ();
if (subFileName != null && !thisEntry.equals (subFileName)) continue;
if (subFileName != null) htParams.put ("subFileName", subFileName);
if (thisEntry.startsWith ("JmolManifest") || haveManifest && exceptFiles == manifest.indexOf ("|" + thisEntry + "|") >= 0) continue;
var bytes = JU.Rdr.getLimitedStreamBytes (zis, ze.getSize ());
if (JU.Rdr.isGzipB (bytes)) bytes = JU.Rdr.getLimitedStreamBytes (zpt.getUnGzippedInputStream (bytes), -1);
if (JU.Rdr.isZipB (bytes) || JU.Rdr.isPngZipB (bytes)) {
var bis = JU.Rdr.getBIS (bytes);
var zipDir2 = zpt.getZipDirectoryAndClose (bis, "JmolManifest");
bis = JU.Rdr.getBIS (bytes);
var atomSetCollections = this.getAtomSetCollectionOrBufferedReaderFromZip (vwr, adapter, bis, fileName + "|" + thisEntry, zipDir2, htParams, ++subFilePtr, asBufferedReader);
if (Clazz.instanceOf (atomSetCollections, String)) {
if (ignoreErrors) continue;
return atomSetCollections;
} else if (Clazz.instanceOf (atomSetCollections, J.adapter.smarter.AtomSetCollection) || Clazz.instanceOf (atomSetCollections, JU.Lst)) {
if (haveManifest && !exceptFiles) htCollections.put (thisEntry, atomSetCollections);
 else vCollections.addLast (atomSetCollections);
} else if (Clazz.instanceOf (atomSetCollections, java.io.BufferedReader)) {
if (doCombine) zis.close ();
return atomSetCollections;
} else {
if (ignoreErrors) continue;
zis.close ();
return "unknown zip reader error";
}} else if (JU.Rdr.isPickleB (bytes)) {
var bis = JU.Rdr.getBIS (bytes);
if (doCombine) zis.close ();
return bis;
} else {
var sData;
if (JU.Rdr.isCompoundDocumentB (bytes)) {
var jd = J.api.Interface.getInterface ("JU.CompoundDocument", vwr, "file");
jd.setStream (zpt, JU.Rdr.getBIS (bytes), true);
sData = jd.getAllDataFiles ("Molecule", "Input").toString ();
} else {
sData = JU.Rdr.fixUTF (bytes);
}var reader = JU.Rdr.getBR (sData);
if (asBufferedReader) {
if (doCombine) zis.close ();
return reader;
}var fname = fileName + "|" + ze.getName ();
ret = adapter.getAtomSetCollectionFromReader (fname, reader, htParams);
if (!(Clazz.instanceOf (ret, J.adapter.smarter.AtomSetCollection))) {
if (ignoreErrors) continue;
zis.close ();
return "" + ret;
}if (haveManifest && !exceptFiles) htCollections.put (thisEntry, ret);
 else vCollections.addLast (ret);
var a = ret;
if (a.errorMessage != null) {
if (ignoreErrors) continue;
zis.close ();
return a.errorMessage;
}}}
if (doCombine) zis.close ();
if (haveManifest && !exceptFiles) {
var list = JU.PT.split (manifest, "|");
for (var i = 0; i < list.length; i++) {
var file = list[i];
if (file.length == 0 || file.indexOf ("#") == 0) continue;
if (htCollections.containsKey (file)) vCollections.addLast (htCollections.get (file));
 else if (JU.Logger.debugging) JU.Logger.debug ("manifested file " + file + " was not found in " + fileName);
}
}if (!doCombine) return vCollections;
var result = (vCollections.size () == 1 && Clazz.instanceOf (vCollections.get (0), J.adapter.smarter.AtomSetCollection) ? vCollections.get (0) :  new J.adapter.smarter.AtomSetCollection ("Array", null, null, vCollections));
if (result.errorMessage != null) {
if (ignoreErrors) return null;
return result.errorMessage;
}if (nFiles == 1) selectedFile = 1;
if (selectedFile > 0 && selectedFile <= vCollections.size ()) return vCollections.get (selectedFile - 1);
return result;
} catch (e$$) {
if (Clazz.exceptionOf (e$$, Exception)) {
var e = e$$;
{
if (ignoreErrors) return null;
JU.Logger.error ("" + e);
return "" + e;
}
} else if (Clazz.exceptionOf (e$$, Error)) {
var er = e$$;
{
JU.Logger.errorEx (null, er);
return "" + er;
}
} else {
throw e$$;
}
}
}, "JV.Viewer,J.api.JmolAdapter,java.io.InputStream,~S,~A,java.util.Map,~N,~B");
Clazz.defineMethod (c$, "clearAndCachePngjFile", 
function (jmb, data) {
jmb.fm.pngjCache =  new java.util.Hashtable ();
if (data == null || data[0] == null) return false;
data[0] = JU.Rdr.getZipRoot (data[0]);
var shortName = J.io.JmolUtil.shortSceneFilename (data[0]);
var cache = jmb.fm.pngjCache;
try {
data[1] = jmb.fm.vwr.getJzt ().cacheZipContents (JU.Rdr.getPngZipStream (jmb.fm.getBufferedInputStreamOrErrorMessageFromName (data[0], null, false, false, null, false, true), true), shortName, cache, false);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
return false;
} else {
throw e;
}
}
if (data[1] == null) return false;
var bytes = data[1].getBytes ();
System.out.println ("jmolutil caching " + bytes.length + " bytes as " + jmb.fm.getCanonicalName (data[0]));
cache.put (jmb.fm.getCanonicalName (data[0]), bytes);
if (shortName.indexOf ("_scene_") >= 0) {
cache.put (J.io.JmolUtil.shortSceneFilename (data[0]), bytes);
bytes = cache.remove (shortName + "|state.spt");
if (bytes != null) cache.put (J.io.JmolUtil.shortSceneFilename (data[0] + "|state.spt"), bytes);
}for (var key, $key = cache.keySet ().iterator (); $key.hasNext () && ((key = $key.next ()) || true);) System.out.println (key);

return true;
}, "J.io.JmolBinary,~A");
Clazz.overrideMethod (c$, "getCachedPngjBytes", 
function (jmb, pathName) {
if (pathName.startsWith ("file:///")) pathName = "file:" + pathName.substring (7);
JU.Logger.info ("JmolUtil checking PNGJ cache for " + pathName);
var shortName = J.io.JmolUtil.shortSceneFilename (pathName);
if (jmb.fm.pngjCache == null && !this.clearAndCachePngjFile (jmb,  Clazz.newArray (-1, [pathName, null]))) return null;
var cache = jmb.fm.pngjCache;
var isMin = (pathName.indexOf (".min.") >= 0);
if (!isMin) {
var cName = jmb.fm.getCanonicalName (JU.Rdr.getZipRoot (pathName));
if (!cache.containsKey (cName) && !this.clearAndCachePngjFile (jmb,  Clazz.newArray (-1, [pathName, null]))) return null;
if (pathName.indexOf ("|") < 0) shortName = cName;
}if (cache.containsKey (shortName)) {
JU.Logger.info ("FileManager using memory cache " + shortName);
return jmb.fm.pngjCache.get (shortName);
}if (!isMin || !this.clearAndCachePngjFile (jmb,  Clazz.newArray (-1, [pathName, null]))) return null;
JU.Logger.info ("FileManager using memory cache " + shortName);
return cache.get (shortName);
}, "J.io.JmolBinary,~S");
Clazz.overrideMethod (c$, "spartanFileList", 
function (zpt, name, type) {
var dirNums = J.io.JmolUtil.getSpartanDirs (type);
if (dirNums.length == 0 && name.endsWith (".spardir.zip") && type.indexOf (".zip|output") >= 0) {
var sname = name.$replace ('\\', '/');
var pt = name.lastIndexOf (".spardir");
pt = sname.lastIndexOf ("/");
sname = name + "|" + name.substring (pt + 1, name.length - 4);
return  Clazz.newArray (-1, ["SpartanSmol", sname, sname + "/output"]);
}return J.io.JmolUtil.getSpartanFileList (name, dirNums);
}, "javajs.api.GenericZipTools,~S,~S");
Clazz.overrideMethod (c$, "getImage", 
function (vwr, fullPathNameOrBytes, echoName) {
var image = null;
var info = null;
var apiPlatform = vwr.apiPlatform;
var createImage = false;
var fullPathName = "" + fullPathNameOrBytes;
if (Clazz.instanceOf (fullPathNameOrBytes, String)) {
var isBMP = fullPathName.toUpperCase ().endsWith ("BMP");
if (fullPathName.indexOf ("|") > 0 || isBMP) {
var ret = vwr.fm.getFileAsBytes (fullPathName, null);
if (!JU.AU.isAB (ret)) return "" + ret;
image = (vwr.isJS ? ret : apiPlatform.createImage (ret));
} else if (vwr.isJS) {
} else if (JU.OC.urlTypeIndex (fullPathName) >= 0) {
try {
image = apiPlatform.createImage ( new java.net.URL (Clazz.castNullAs ("java.net.URL"), fullPathName, null));
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
return "bad URL: " + fullPathName;
} else {
throw e;
}
}
} else {
createImage = true;
}} else if (vwr.isJS) {
image = fullPathNameOrBytes;
} else {
createImage = true;
}if (createImage) image = apiPlatform.createImage ("\1close".equals (fullPathNameOrBytes) ? "\1close" + echoName : fullPathNameOrBytes);
{
info = [echoName, fullPathNameOrBytes];
}try {
if (!apiPlatform.waitForDisplay (info, image)) return null;
{
return null;
}} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
return e.toString () + " opening " + fullPathName;
} else {
throw e;
}
}
}, "JV.Viewer,~O,~S");
});

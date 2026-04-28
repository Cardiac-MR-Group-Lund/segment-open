function tf = myisdicom(filename)
% function to check if a file is a DICOM, but excludes files with certain
% extensions as non-dicom right away
% files no on this list are checked with Matlab own isdicom

% get file part and only process files that are not on exclusion list
[~,~,ext] = fileparts(filename);
isnotadicom = shouldignorefile(ext);

if isnotadicom
  % return false if file extension is on exclusion list
  tf = false;
else
  % check this file with Matlab own isdicom
  tf = isdicom(filename);
end


function tf = shouldignorefile(ext)
% load exclusion list
persistent excllist addonlist
if isempty(excllist)
  excllist = getexclusionlist;
end

if isempty(addonlist)
  addonlist = getexcladdonlist;
end
tf = any(strcmpi(ext,excllist),'all') || any(contains(ext,addonlist,'IgnoreCase',true));

function addonlist = getexcladdonlist
% this is a list of extensions that can be followed by digits
addonlist = {
  '.data';
  '.dll';  
  '.exe';
  '.fon';
  '.hBaked';
  '.ini';
  '.lic';
  '.log';
  '.mex';
  '.mui';
  '.ppkg';
  '.ps';
  '.ttf';
  '.sys';
  '.vc';
  
  };
function excllist = getexclusionlist

excllist = {
  
'.3g2'; ... 3GPP2 multimedia file
'.3gp'; ...3GPP multimedia file
'.7z'; ... 7-Zip compressed file

'.acl'; ...AutoCorrect List File
'.afm'; 
'.aft';
'.ai'; ... Adobe Illustrator file
'.AIFF' ;'.AIF' ;	...Audio Interchange File Format
'.alco';...
'.amp'; ...
'.ansi'; ...Mistaken TXT Document Format
'.apk'; ... Android package file
'.appx';...
'.appxbundle';...
'.arj'; ... ARJ compressed file
'.asc'; ...armored ASCII file
'.ascii'; ... ASCII Text Format
'.asd'; ... Microsoft Word Auto recovery Document
'.asp'; '.aspx'; ... Active Server Page file
'.asv'; ... Matlabs backup
'.AU'; 	...Basic Audio
'.avastlic'; ...
'.AVI'; ...	Multimedia Audio/Video

'.bak'; ... Backup file
'.BAT'; ... PC batch file
'.bin'; ... Binary disc image
'.blf'; ...
'.BMP'; 	...Windows BitMap
'.bto'; ...
'.bytecode'; ...

'.c';    ... cfiles
'.cab'; ... Windows Cabinet file
'.cat'; ...
'.cda';  ...CD audio track file
'.cer'; ... Internet security certificate
'.cfg'; ... Configuration file
'.cfm'; ... ColdFusion Markup file
'.cfx'; ... flowcal file
'.cgi';  '.pl'; ... Perl script file
'.CHK'; ... files are used by ScanDisk, a disk diagnostics and repair tool.
'.cookie'; ...
'.com'; ... MS-DOS command file
'.conf'; ... configuration file
'.config'; ... configuration file
'.CLASS' ; '.JAVA';	...Java files
'.clx' ; 	... file used by Adobe
'.cmd'; ...
'.cpa'; ...
'.cpl'; ... Windows Control panel file
'.cpp'; ... C++ source code file
'.cs'; ... Visual C# source code file
'.csf'; ... Color management settings file used by Adobe
'.css'; ... Cascading Style Sheet file
'.CSV'; 	...Comma separated, variable length file (Open in Excel)
'.customDestinations-ms'; ...
'.cur'; ... Windows cursor file
'.cva'; ... 
'.CVS'; 	...Canvas

'.dat'; ... Data file
'.data'; ...DATA files mostly belong to Analysis Studio
'.db' ; '.dbf'; '.dbt'; ... Database file
'.deb'; ... Debian software package file
'.DIF'; 	...Data Interchange format
'.dmg'; ... macOS X disk image
'.dmp'; ... Dump file
'.dll';  ... dynamic library
'.DOC'; '.DOCX'; 	...Microsoft Word for Windows/Word97
'.drv'; ... Device driver file
'.dsp'; ...
'.dt'; ... Data Files used in operating system Windows

'.edb'; ... Database created by Microsoft Exchange Server
'.ehsc'; ...
'.ens'; ...
'.email'; ... Outlook Express e-mail message file.
'.eml'; ... E-mail message file from multiple e-mail clients, including Gmail.
'.EPS';	...Encapsulated PostScript
'.err'; ... error log files
'.etl'; ...Microsoft Tracelog
'.eth'; ...
'.evtx'; ...  Vista Event Log file
'.EXE'; 	...PC Application

'.fig';  ... Matlab figure file
'.flv';   ... Adobe Flash file
'.FM3'; 	...Filemaker Pro databases (the numbers following represent the version #)
'.fnt'; ...Windows font file
'.fon'; ... Generic font file
'.ft'; ...

'.gadget'; ... Windows gadget
'.geo'; ...
'.GIF'; 	...Graphics Interchange Format
'.gz'; ... Tarball compressed file

'.h'; ... C, C++, and Objective-C header file
'.h264';  ... H.264 video file
'.hds'; ... Parallels Desktop virtualization file
'.hpps'; ...
'.hsc'; ...
'.HQX'; 	...Macintosh BinHex
'.htf'; ...
'.HTM'; '.HTML'; 	...Web page source text
'.hve'; ...Windows Registry Hive File

'.icc'; ... color profile format standardized by the International Color Consortium (ICC)
'.icns'; ...  macOS X icon resource file
'.ico'; ... Icon file
'.index';...
'.inf'; ...
'.ini'; ... Initialization file
'.iso'; ... ISO disc image

'.jar'; ... Java Archive file
'.jfm'; ... JetForm Filler, a software for filling on-screen forms
'.job'; ...
'.JPG' ; 'JPEG'; 	...JPEG graphic
'.jrs'; ... Transaction log file created by Microsoft Exchange
'.js'; ... JavaScript file
'.json'; ... JSON files
'.jsp'; ... Java Server Page file
'.jtx';

'.key';  ...Keynote presentation
'.kic'; ...

'.ldb'; ... level db
'.lib';  ...Library 
'.lic';  ...Medvisos file
'.lnk'; ... Windows shortcut file
'.lock'; ...
'.log'; ...  log file
'.log1'; '.log2'; ...  log file
'.lproj'; ... Localized Project Folder
'.lut'; ... look up table

'.m'; ... Matlab function files
'.m3u';
'.m4a'; 
'.m4v'; ... Apple MP4 video file
'.MAC'; 	...MacPaint
'.MAP'; 	...Web page imagemap
'.mat'; ... Matlab file
'.MDB'; 	...MS Access database
'.mexa64'; '.mexmaci64'; '.mexw64'; '.mexw32'; ...Matlabs mex files
'.mf';
'.MID' ; '.MIDI'; 	...MIDI sound
'.mine'; ... SVN file
'.mkv'; ... Matroska Multimedia Container
'.mof'; ... Managed Object Format
'.MOV' ; '.QT'; 	...QuickTime Audio/Video
'.msi'; ... Windows installer package
'.msg'; ... Microsoft Outlook e-mail message file.
'.mss'; ...Microprocessor Software Specification File
'.mp3';  ...MP3 audio file
'.mp4' ; ... MPEG4 video file
'.mpa';  ...  MPEG-2 audio file
'.mpg' ; '.mpeg'; ... MPEG video file
'.MTB';  '.MTW'; 	...MiniTab
'.mui'; ...Multilingual User Interface component

'.node'; ...
'.nrm'; ...
'.nv'; ... Medvisos file
'.nvph';...

'.obj'; ... Geometry of 3D objects
'.odl'; ... File written in the Object Description Language
'.odp'; ... OpenOffice Impress presentation file
'.ods'; ... OpenOffice Calc spreadsheet file
'.odt'; ... OpenOffice Writer document file
'.oec'; ...
'.oft'; ... Microsoft Outlook e-mail template file.
'.ogg'; ... Ogg Vorbis audio file
'.old'; ...
'.ORT';...  Oric Raw Tape format file.
'.ost'; ... Microsoft Outlook offline e-mail storage file.
'.otf'; ... Open type font file

'.p';    ...Matlab protected file
'.P65';  ...PageMaker (the numbers following represent the version #) P=publication, T=template
'.p7s'; ...
'.pak'; ... an archive used by video games
'.par'; ... medvisos par file
'.part'; ... Partially downloaded file
'.patch'; ...File created by Mercurial, a control management developer tool;
'.pcb'; ... associated with the Microsoft PowerPoint presentation 
'.PDF';  ...Acrobat -Portable document format
'.pfb'; ...
'.pfm'; ...
'.phm';
'.php'; ... PHP file
'.pkg'; ... Package file
'.PNG'; ...Portable Network Graphics
'.pol'; ...
'.PPT';  '.PPTX'; 	...PowerPoint
'.pps'; ... PowerPoint slide show
'.prev'; ...
'.pri'; ...
'.provxml';...
'.ps'; ... PostScript file
'.PSD'; 	...Adobe PhotoShop
'.PSP'; 	...PaintShop Pro
'.pst'; ... Microsoft Outlook e-mail storage file.
'.py'; ... Python file

'.QXD'; 	...QuarkXPress

'.RA'; 	...RealAudio
'.rar'; ... RAR file
'.rom';...
'.rm'; ... RealMedia file
'.rpm'; ... Red Hat Package Manager
'.rss'; ... RSS file
'.RTF'; 	...Rich Text Format

'.sav'; ... Save file (e.g., game save file)
'.sdb'; ... Windows registry database files
'.sh'; ... Bash shell script
'.SIT'; 	...Stuffit Compressed Archive
'.sql'; ... SQL database file
'.STL'; ... 	Stereolithography
'.svg'; ... Scalable Vector Graphics file
'.swf'; ... Shockwave flash file
'.swift'; ... Swift source code file
'.svn-base'; ... SVN file
'.sys'; ... Windows system file

'.T65'; 	...PageMaker (the numbers following represent the version #) P=publication, T=template
'.table'; ...
'.TAR'; 	...UNIX TAR Compressed Archive
'.tex';  ...  A LaTeX document file
'.tf'; ... 
'.TFM'; ... FormTool Gold form file
'.th'; ...
'.TIF'; 	...TIFF graphic
'.tmp'; ... Temporary file
'.toast'; ... Toast disc image
'.toc'; ... associated with the Eudora email
'.ttf'; ... TrueType font file
'.tv13';
'.TXT'; 	...ASCII text (Mac text does not contain line feeds--use DOS Washer Utility to fix)

'.vb'; ... Visual Basic file
'.vcd'; ... Virtual CD
'.vcf'; ... E-mail contact file.
'.vdm'; ... Windows defender definition updates
'.vf'; ...
'.vob'; ... DVD Video Object
'.vol'; ...
'.vp'; ...
'.vtk'; ...  Visualization Toolkit format file

'.WAV'; 	...Windows sound
'.WK3'; 	...Lotus 1-2-3 (the numbers following represent the version #)
'.WKS'; 	...MS Works
'.wma'; ...WMA audio file
'.wmv'; ... Windows Media Video file
'.WPD';  '.WP5'; ...	WordPerfect (the numbers following represent the version #)
'.wpl'; ... Windows Media Player playlist
'.wsf'; ... Windows Script File

'.xhtml'; ... XHTML file
'.XLS';  '.XLSX'; ... 	Excel spreadsheet
'.xlsb'; ... Microsoft Excel binary worksheet
'.xlsm'; ... Microsoft Excel file with macros
'.xml'; ... XML file
'.xml1'; ...
'.xsd'; ... XML Schema Description

'.ytr'; ... an IRIS OCR data file

'.z'; ... Z compressed file
'.ZIP'... 	PC Zip Compressed Archive 

};
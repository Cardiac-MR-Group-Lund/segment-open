function segloaderprogressbar(action, event)
  
  persistent wb;
  
  switch action
    case 'init'
      wb = segloaderprogressbarinit();
    case 'close'
      if not(isempty(wb))
        segloaderprogressbarclose(wb);
        wb = [];
      end
    case 'update'
      if not(isempty(wb))
        wb = segloaderprogressbarupdate(wb, event);
      end
  end
end

function wb = segloaderprogressbarinit()
  wb = struct();
  wb.h = waitbar(0,dprintf('Loading files...'));
  wb.filesToRead = [];
  wb.filesRead = 0;
  wb.lastWbValue = 0;
  wb.numStacks = [];
  wb.stacksRendered = 0;
  wb.numDicomsUniqueLines = [];
  wb.numDicomsDoneUniqueLines = 0;
  
  wb.numDicomsRmDups = [];
  wb.numDicomsDoneRmDups = 0;
  
  wb.numDicomsSetDcm = [];
  wb.numDicomsDoneSetDcm = 0;

  wb.numImagesImData = [];
  wb.numImagesDoneImData = 0;
end

function segloaderprogressbarclose(wb)
  close(wb.h);
end

function wb = segloaderprogressbarupdate(wb, event)
  READFILESSIZE = 20;
  HASSEGMENTDATA1SIZE = 2;
  HASSEGMENTDATA2SIZE = 1;
  IGNOREMESIZE = 2;
  UNIQUELINESSIZE = 5;
  RENDERSIZE = 70;
  
  REMOVEDUPSSIZE = 5/100;
  SETDICOMSSIZE = 25/100;
  MAKEIMBASESIZE = 5/100;
  MAKEIMDATASIZE = 25/100;
  RAWSTACK1SIZE = 5/100;
  RAWSTACK2SIZE = 5/100;
  RAWSTACK3SIZE = 5/100;
  RAWSTACK4SIZE = 5/100;
  RAWSTACK5SIZE = 5/100;
  RAWSTACK6SIZE = 5/100;
  RAWSTACK7SIZE = 5/100;
  RAWSTACK8SIZE = 5/100;
  
  switch event.name
    case 'readfiles'
      wb.filesToRead = event.numfiles;
    case 'readfile'
      wb.filesRead = wb.filesRead + 1;
      wb = segloaderprogressbarupdatewb( ...
        wb.filesRead / wb.filesToRead * READFILESSIZE, wb);
    case 'readmat'
      wb = segloaderprogressbarupdatewb(READFILESSIZE, wb);
    case 'rendermat'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+RENDERSIZE, wb);
    case 'hassegmentdata1'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE, wb);
    case 'hassegmentdata2'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE, wb);
    case 'ignoreme'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE, wb);
    case 'uniquelinesstart'
      wb.numDicomsUniqueLines = event.numdicoms;
      wb.numDicomsDoneUniqueLines = 0;
    case 'uniqueline'
      wb.numDicomsDoneUniqueLines = wb.numDicomsDoneUniqueLines + 1;
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE + ...
        UNIQUELINESSIZE * ...
        (wb.numDicomsDoneUniqueLines/wb.numDicomsUniqueLines), wb);
    case 'renderstacksstart'
      wb.numStacks = event.numstacks;
      wb.stacksRendered = 0;
    case 'renderstack'
      wb.stacksRendered = wb.stacksRendered + 1;
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered) / wb.numStacks , ...
        wb);
    case 'removedupsstart'
      wb.numDicomsRmDups = event.numdicoms;
      wb.numDicomsDoneRmDups = 0;
    case 'removedup'
      wb.numDicomsDoneRmDups = wb.numDicomsDoneRmDups + 1;
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + (wb.numDicomsDoneRmDups/wb.numDicomsRmDups)*REMOVEDUPSSIZE) / wb.numStacks , ...
        wb);
    case 'setdicomsstart'
      wb.numDicomsSetDcm = event.numdicoms;
      wb.numDicomsDoneSetDcm = 0;
    case 'setdicom'
      wb.numDicomsDoneSetDcm = wb.numDicomsDoneSetDcm + 1;
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + (REMOVEDUPSSIZE + (wb.numDicomsDoneSetDcm/wb.numDicomsSetDcm)*SETDICOMSSIZE)) / wb.numStacks , ...
        wb);
    case 'makeimbase'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + (REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE)) / wb.numStacks , ...
        wb);
    case 'makeimdatastart'
      wb.numImagesImData = event.numimages;
      wb.numImagesDoneImData = 0;
    case 'makeimdata'
      wb.numImagesDoneImData = wb.numImagesDoneImData + 1;
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + (REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + (wb.numImagesDoneImData/wb.numImagesImData)*MAKEIMDATASIZE)) / wb.numStacks , ...
        wb);
    case 'rawstack1'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack2'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack3'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE + RAWSTACK3SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack4'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE + RAWSTACK3SIZE + RAWSTACK4SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack5'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE + RAWSTACK3SIZE + RAWSTACK4SIZE + RAWSTACK5SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack6'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE + RAWSTACK3SIZE + RAWSTACK4SIZE + RAWSTACK5SIZE + RAWSTACK6SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack7'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE + RAWSTACK3SIZE + RAWSTACK4SIZE + RAWSTACK5SIZE + RAWSTACK6SIZE + RAWSTACK7SIZE) / wb.numStacks , ...
        wb);
    case 'rawstack8'
      wb = segloaderprogressbarupdatewb(READFILESSIZE+HASSEGMENTDATA1SIZE+HASSEGMENTDATA2SIZE+IGNOREMESIZE+UNIQUELINESSIZE + ...
        RENDERSIZE*(wb.stacksRendered + REMOVEDUPSSIZE + SETDICOMSSIZE + MAKEIMBASESIZE + MAKEIMDATASIZE + RAWSTACK1SIZE + RAWSTACK2SIZE + RAWSTACK3SIZE + RAWSTACK4SIZE + RAWSTACK5SIZE + RAWSTACK6SIZE + RAWSTACK7SIZE + RAWSTACK8SIZE) / wb.numStacks , ...
        wb);
  end
end

function wb = segloaderprogressbarupdatewb(v, wb)
  v = floor(v);
  if v == wb.lastWbValue
    return
  end
  try
    waitbar(v / 100, wb.h);
  catch
    %ignore
  end
  wb.lastWbValue = v;
end
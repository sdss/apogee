pro apql_listensocket,astring,file=file

if n_elements(file) eq 0 then file='apql_messages.txt'

; Read the message text file and check for
; a new message
READLINE,file,lines,count=nlines,comment='#'
if nlines eq 0 then begin
  print,'No message'
  return
endif

; Create system variable if it doesn't exist
DEFSYSV,'!message',exists=message_exists
if not message_exists then begin
  message = {nlines:-1L,messages:PTR_NEW(),lastcount:0L}
  message.nlines = nlines
  message.messages = PTR_NEW(lines)
  message.lastcount = -1
  DEFSYSV,'!message',message
endif

; Add new messages
!message.nlines = nlines
!message.messages = PTR_NEW(lines)

; Get next message
count = !message.lastcount+1
if count gt !message.nlines-1 then return ; no new messages
astring = (*!message.messages)[count]

; Update lastcount
!message.lastcount++

;stop

end

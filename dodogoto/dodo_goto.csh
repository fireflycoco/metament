#!/bin/csh -f

set dest = `dodo_goto_parser ${*}`

if ( "${dest}" != "." ) then
	cd "$dest"
endif

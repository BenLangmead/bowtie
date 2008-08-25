 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: index_esa_drawing.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_ESA_DRAWING_H
#define SEQAN_HEADER_INDEX_ESA_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TFile, typename TText, typename TESASpec>
void write(TFile & file, 
	   Index<TText, Index_ESA<TESASpec> > & stree,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Index<TText, Index_ESA<TESASpec> > TIndex;
	
	_streamWrite(file, "digraph G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TIndex, TopDown<ParentLinks<Preorder> > >::Type TIterator;
	typedef typename Iterator<TIndex, TopDown<> >::Type TIteratorSimple;
	TIterator it(stree);

	for(;!atEnd(it);++it) 
	{
		// dump node
       		_streamWrite(file, "\"[");
 		_streamPutInt(file, value(it).range.i1);
		_streamPut(file, ':');
		_streamPutInt(file, value(it).range.i2);
       		_streamWrite(file, ")\"");
       		if (!isRightTerminal(it))
			_streamWrite(file, " [style = dashed]");
       		_streamWrite(file, ";\n");

		// dump edge from parent (if not root)
		if (!isRoot(it)) {
			TIteratorSimple src(container(it), nodeUp(it));

			_streamWrite(file, "\"[");
			_streamPutInt(file, value(src).range.i1);
			_streamPut(file, ':');
			_streamPutInt(file, value(src).range.i2);
			_streamWrite(file, ")\"");

			_streamWrite(file, " -> ");

			_streamWrite(file, "\"[");
			_streamPutInt(file, value(it).range.i1);
			_streamPut(file, ':');
			_streamPutInt(file, value(it).range.i2);
			_streamWrite(file, ")\"");

			_streamWrite(file, " [label = \"");
			_streamWrite(file, parentEdgeLabel(it));
			_streamWrite(file, "\"];\n");
		}
	}
	_streamPut(file, '\n');

	_streamWrite(file, "}\n");
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

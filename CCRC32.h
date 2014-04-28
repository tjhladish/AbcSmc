///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Copyright © NetworkDLS 2002, All rights reserved
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF 
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A 
// PARTICULAR PURPOSE.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _CCRC32_H
#define _CCRC32_H
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CCRC32{

	public:
		CCRC32(void);
		~CCRC32(void);

		void Initialize(void);

		bool FileCRC(const char *sFileName, unsigned int *iOutCRC);
		bool FileCRC(const char *sFileName, unsigned int *iOutCRC, size_t iBufferSize);
		unsigned int FileCRC(const char *sFileName);
		unsigned int FileCRC(const char *sFileName, size_t iBufferSize);

		unsigned int FullCRC(const unsigned char *sData, size_t iDataLength);
		void FullCRC(const unsigned char *sData, size_t iLength, unsigned int *iOutCRC);

		void PartialCRC(unsigned int *iCRC, const unsigned char *sData, size_t iDataLength);

	private:
		unsigned int Reflect(unsigned int iReflect, const char cChar);
		unsigned int iTable[256]; // CRC lookup table array.
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

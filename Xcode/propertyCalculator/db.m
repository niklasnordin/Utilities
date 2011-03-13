//
//  db.m
//  propertyCalculator
//
//  Created by Niklas Nordin on 2011-01-03.
//  Copyright 2011 nequam. All rights reserved.
//

#import "db.h"
#import "specie.h"

@implementation db


-(void) create
{
	// set root_ and append the trailing '/' if it is missing
	root_ = @"/Users/niklasnordin/Documents/CoreData/propertyCalculator/Docs/";
	if (![root_ hasSuffix:@"/"])
	{
		root_ = [root_ stringByAppendingFormat:@"/"];
	}
	
	files_ = [[NSMutableArray alloc] init];
	species_ = [[NSMutableArray alloc] init];
	
	//NSLog(@"files_ count = %d\n", [files_ count]);
	NSFileManager *manager_ = [NSFileManager defaultManager];
	
	// get all the files in the root_ directory with the correct '.prop' suffix
	NSArray *children = [manager_ contentsOfDirectoryAtPath: root_ error:NULL];
	for (NSString *filename in children)
	{
		if ([filename hasSuffix: @".prop"])
		{
			[files_ addObject: [root_ stringByAppendingFormat: @"%@", filename] ];
		}
	}
	
	//NSLog(@"files_ count = %d\n", [files_ count]);
	// read all the .prop files 
	for (NSString *filename in files_)
	{
		NSString *fileString = [NSString stringWithContentsOfFile: [filename stringByStandardizingPath] encoding:NSUTF8StringEncoding error:NULL];
		NSString *name = [self getSpecieName: fileString];
		
		// check if the species name is already in there...get the object.
		int size = [species_ count];
		NSLog(@"size = %d.",size);
		
		if (size > 0)
		{
			// check if species already exist
			int index = -1;
			for(int i=0; i<size; i++)
			{
				NSString *n = [[species_ objectAtIndex:i] name];
				if ([name isEqualToString: n])
				{
					index = i;
				}
			}
			
			// specie already exists
			if (index < 0)
			{
				NSLog(@"*** adding specie %@",name);
				specie *s = [[specie alloc] init];
				[s setName: name];
				[species_ addObject: s];
				NSArray *k = [self getProperties: fileString];
				for(NSString *a in k)
				{
					NSLog(@"--------");
					NSLog(@"%@", a);
				}
				NSLog(@"*** added specie %@",name);
			}
			else 
			{
				NSLog(@"<<< specie already exists %@",name);
			}

			NSLog(@"index = %d", index);
		}
		else 
		{
			// no available species, add the first one
			NSLog(@"adding specie %@",name);
			specie *s = [[specie alloc] init];
			[s setName: name];
			[species_ addObject: s];
			
			// and allocate the property list
			NSMutableArray *properties = [s properties];
			properties = [[NSMutableArray alloc] init];

			NSLog(@"added specie %@",name);
			NSArray *propsFromFile = [self getProperties: fileString];
			int nProps = [propsFromFile count];
			for(int i=1; i<nProps; i++)
			//for(NSString *a in props)
			{
				NSString *a = [propsFromFile objectAtIndex:i];
				NSLog(@"--------");
				NSLog(@"%@", a);
				// check previous properties
				NSMutableArray *prevProps = [s properties];
				int nPrev = [prevProps count];
				NSLog(@"nPrev=%d", nPrev);
				[properties addObject: a];
			}
		}
	}
	
}

-(NSMutableArray *) species
{
	return species_;
}

-(NSString *) getSpecieName:(NSString *)fileString
{
	NSArray *lines = [fileString componentsSeparatedByString:@"\n"];
	int nLines = [lines count];
	bool found = FALSE;
	int i=0;
	NSString *name = @"";
	while((i<nLines) && (!found))
	{
		NSArray *words = [[lines objectAtIndex: i] componentsSeparatedByString:@" "];
		NSString *identifier = [words objectAtIndex:0];
		//NSLog(@"%@", identifier);
		if ([identifier isEqualToString: @"name"])
		{
			name = [words objectAtIndex:1];
			found = TRUE;
		}
		i++;
	}
	return name;
}

-(NSArray *) getProperties:(NSString *)fileString
{
	NSArray *props = [fileString componentsSeparatedByString: @"property "];
	//int nProps = [props count];
	//NSLog(@"nProps = %d.", nProps);
	return props;
}

@end

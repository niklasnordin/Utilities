//
//  db.h
//  propertyCalculator
//
//  Created by Niklas Nordin on 2011-01-03.
//  Copyright 2011 nequam. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface db
: NSObject 
{

	NSString *root_;
	NSMutableArray *files_;
	NSMutableArray *species_;
	
}

-(void) create;

-(NSMutableArray *) species;

-(NSString *) getSpecieName:(NSString *) fileString;
-(NSArray *) getProperties:(NSString *) fileString;

@end

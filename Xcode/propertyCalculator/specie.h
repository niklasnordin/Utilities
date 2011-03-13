//
//  specie.h
//  propertyCalculator
//
//  Created by Niklas Nordin on 2011-01-05.
//  Copyright 2011 nequam. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface specie
: NSObject
{
	NSString *name_;
	NSMutableArray *properties_;
	NSMutableArray *function_;
	NSMutableArray *coefficients_;
}

//-(void) init;

-(NSString *) name;
-(void) setName:(NSString *) name;

-(NSMutableArray *) properties;

@end
